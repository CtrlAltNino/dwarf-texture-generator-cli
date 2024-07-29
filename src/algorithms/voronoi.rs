use clap::Args;
extern crate hilbert;
extern crate nalgebra as na;
extern crate rand;
use delaunator;

use crate::algorithms::Noise;

use crate::rector::*;
use core::{f64, panic};
use hilbert::point as hpoint;
use na::{ComplexField, Const, OPoint, Point, Point2, RealField, Vector2};
use rand::Rng;
use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::fmt;
use std::hash::{Hash, Hasher};
use std::io::Write;
use std::ops::{Add, Div, Mul, Sub};
use std::time::Instant;
use std::{
    borrow::Borrow,
    cell::RefCell,
    collections::HashSet,
    rc::{Rc, Weak},
};

#[derive(Args, Debug)]
pub struct Voronoi {
    // #[arg(short, long)]
    // seed: u32,
    #[arg(short, long)]
    cells: u32,
}

extern crate image;
#[derive(Debug, Clone)]
struct Edge {
    start: Point2<f64>,
    end: Point2<f64>,
}

#[derive(Debug, Clone)]
struct Triangle {
    vertices: [Point2<f64>; 3],
    edges: [Edge; 3],
    circumcenter: Point2<f64>,
    neighbours: RefCell<Vec<Weak<Triangle>>>,
    checked: RefCell<bool>,
}

fn smallest_angle(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> f64 {
    // Helper function to calculate the distance between two points
    fn distance(p1: (f64, f64), p2: (f64, f64)) -> f64 {
        ((p2.0 - p1.0).powi(2) + (p2.1 - p1.1).powi(2)).sqrt()
    }

    // Calculate the lengths of the sides of the triangle
    let ab = distance(a, b);
    let bc = distance(b, c);
    let ca = distance(c, a);

    // Helper function to calculate the angle opposite to the given side
    fn angle(opposite: f64, side1: f64, side2: f64) -> f64 {
        ((side1.powi(2) + side2.powi(2) - opposite.powi(2)) / (2.0 * side1 * side2)).acos()
    }

    // Calculate all three angles in the triangle
    let angle_a = angle(bc, ab, ca);
    let angle_b = angle(ca, bc, ab);
    let angle_c = angle(ab, ca, bc);

    // Return the smallest angle
    angle_a.min(angle_b).min(angle_c)
}

impl Triangle {
    fn new(p1: Point2<f64>, p2: Point2<f64>, p3: Point2<f64>) -> Self {
        let edges = [
            Edge { start: p1, end: p2 },
            Edge { start: p2, end: p3 },
            Edge { start: p3, end: p1 },
        ];

        let d = 2.0 * (p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y));
        let ux = ((p1.x * p1.x + p1.y * p1.y) * (p2.y - p3.y)
            + (p2.x * p2.x + p2.y * p2.y) * (p3.y - p1.y)
            + (p3.x * p3.x + p3.y * p3.y) * (p1.y - p2.y))
            / d;
        let uy = ((p1.x * p1.x + p1.y * p1.y) * (p3.x - p2.x)
            + (p2.x * p2.x + p2.y * p2.y) * (p1.x - p3.x)
            + (p3.x * p3.x + p3.y * p3.y) * (p2.x - p1.x))
            / d;

        Triangle {
            vertices: [p1, p2, p3],
            edges,
            circumcenter: Point2::new(ux, uy),
            neighbours: RefCell::new(vec![Weak::new(), Weak::new(), Weak::new()]),
            checked: RefCell::new(false),
        }
    }

    fn update_neighbourhood(self: &Rc<Self>, other: &Rc<Self>) -> bool {
        let Triangle { edges: edges1, .. } = Rc::as_ref(self);
        let Triangle { edges: edges2, .. } = Rc::as_ref(other);

        for (i, e) in edges1.iter().enumerate() {
            for (j, f) in edges2.iter().enumerate() {
                if e == f {
                    self.neighbours.borrow_mut()[i] = Rc::downgrade(other);
                    other.neighbours.borrow_mut()[j] = Rc::downgrade(self);
                    return true;
                }
            }
        }
        false
    }

    fn say_goodbye(&self) {
        // Iterate over all neighbours
        for weak_neighbour in self.neighbours.borrow().iter() {
            if let Some(neighbour) = weak_neighbour.upgrade() {
                // Get a mutable reference to the neighbour's neighbours
                let mut neighbours_list = neighbour.neighbours.borrow_mut();

                // Find and remove this triangle from the neighbour's neighbours list
                for neighbour_weak_ref in neighbours_list.iter_mut() {
                    if neighbour_weak_ref
                        .upgrade()
                        .map_or(false, |tri| *tri == *self)
                    {
                        *neighbour_weak_ref = Weak::new();
                    }
                }
            }
        }
    }

    fn contains_point(&self, point: &Point2<f64>) -> bool {
        let d1 = sign(point, &self.vertices[0], &self.vertices[1]);
        let d2 = sign(point, &self.vertices[1], &self.vertices[2]);
        let d3 = sign(point, &self.vertices[2], &self.vertices[0]);

        let has_neg = (d1 < 0.0) || (d2 < 0.0) || (d3 < 0.0);
        let has_pos = (d1 > 0.0) || (d2 > 0.0) || (d3 > 0.0);

        !(has_neg && has_pos)
    }

    fn circumcircle_contains(&self, point: &Point2<f64>) -> bool {
        let a = self.vertices[0];
        let b = self.vertices[1];
        let c = self.vertices[2];

        // Creating the matrix
        let matrix = na::Matrix4::new(
            a.x,
            a.y,
            a.coords.norm_squared(),
            1.0,
            b.x,
            b.y,
            b.coords.norm_squared(),
            1.0,
            c.x,
            c.y,
            c.coords.norm_squared(),
            1.0,
            point.x,
            point.y,
            point.coords.norm_squared(),
            1.0,
        );

        // Calculate the determinant
        let det = matrix.determinant();

        // If the determinant is positive, the point is inside the circumcircle
        det > 0.0
    }
}

impl Hash for Triangle {
    fn hash<H: Hasher>(&self, state: &mut H) {
        let [p1, p2, p3] = self.vertices;
        Hash::hash(
            &[
                p1.x as u64,
                p1.y as u64,
                p2.x as u64,
                p2.y as u64,
                p3.x as u64,
                p3.y as u64,
            ],
            state,
        );
    }
}

impl fmt::Display for Edge {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}<->{})", self.start, self.end)
    }
}

impl fmt::Display for Triangle {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let Triangle {
            edges:
                [Edge { start: p1, end: p2 }, Edge { start: p3, end: p4 }, Edge { start: p5, end: p6 }],
            ..
        } = self;
        write!(f, "({}<->{}<->{})", p1, p2, p4)
    }
}

fn sign(p1: &Point2<f64>, p2: &Point2<f64>, p3: &Point2<f64>) -> f64 {
    (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y)
}

impl PartialEq for Edge {
    fn eq(&self, other: &Self) -> bool {
        (self.start == other.start && self.end == other.end)
            || (self.start == other.end && self.end == other.start)
    }
}

impl Eq for Triangle {}

impl PartialOrd for Triangle {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.vertices.partial_cmp(&other.vertices)
    }
}

impl Ord for Triangle {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        let Triangle {
            vertices: [p1, p2, p3],
            ..
        } = self;
        let Triangle {
            vertices: [q1, q2, q3],
            ..
        } = other;

        fn my_cmp(p: &Point2<f64>, q: &Point2<f64>) -> std::cmp::Ordering {
            match f64::total_cmp(&p.x, &q.x) {
                r @ Ordering::Less | r @ Ordering::Greater => r,
                Ordering::Equal => f64::total_cmp(&p.y, &q.y),
            }
        }

        match my_cmp(&p1, &q1) {
            r @ Ordering::Less | r @ Ordering::Greater => r,
            Ordering::Equal => match my_cmp(&p2, &q2) {
                r @ Ordering::Less | r @ Ordering::Greater => r,
                Ordering::Equal => my_cmp(&p2, &q2),
            },
        }
    }
}

impl PartialEq for Triangle {
    fn eq(&self, other: &Self) -> bool {
        self.vertices.iter().all(|v| other.vertices.contains(v))
    }
}

fn add_neighbours<'a, I>(triangles: I, new_triangle: &Rc<Triangle>)
where
    I: IntoIterator<Item = &'a Rc<Triangle>>,
{
    for t in triangles {
        new_triangle.update_neighbourhood(&t);
    }
}

fn print_progress(p: usize, q: usize) {
    let a = 50;
    let b = (p * a) / q;
    let s = "█".repeat(b);
    let t = "░".repeat(a - b);

    eprint!("\r{}{} {}/{}", s, t, p, q);
    std::io::stderr().flush().unwrap();
}

fn get_region(start: &Rc<Triangle>, point: &Point2<f64>, bad_triangles: &mut Vec<Rc<Triangle>>) {
    // let mut res = vec![Rc::clone(start)];
    bad_triangles.push(Rc::clone(start));

    fn inner(si: &Rc<Triangle>, p: &Point2<f64>, r: &mut Vec<Rc<Triangle>>) {
        si.checked.replace(true);
        for neighbour_weak_ref in si.neighbours.borrow().iter() {
            if let Some(n_ref) = neighbour_weak_ref.upgrade() {
                if !*n_ref.checked.borrow() && n_ref.circumcircle_contains(p) {
                    r.push(Rc::clone(&n_ref));
                    inner(&n_ref, p, r);
                }
            }
        }
    }
    inner(start, point, bad_triangles);

    for r in bad_triangles {
        r.checked.replace(false);
    }
}

fn get_relation(t1: &Rc<Triangle>, t2: &Rc<Triangle>) -> Option<(u8, u8)> {
    for (i, e) in t1.edges.iter().enumerate() {
        for (j, f) in t2.edges.iter().enumerate() {
            if e == f {
                return Some((i as u8, j as u8));
            }
        }
    }
    None
}

fn get_polygon_edge(triangle: &Rc<Triangle>, triangles: &Vec<Rc<Triangle>>) -> Option<u8> {
    let mut relations = vec![];
    for t in triangles {
        if let Some(tup) = get_relation(triangle, t) {
            relations.push(tup);
        }
    }

    if relations.len() == 2 {
        if let [(i, _), (j, _)] = relations[..] {
            let v: Vec<u8> = (0..3).filter(|t| *t != i && *t != j).collect();
            return v.first().copied();
        }
    }
    None
}

#[derive(Debug)]
struct Triangles {
    triangles: Vec<Triangle>,
}

fn scanline_fill<F>(vertices: &[delaunator::Point], mut draw_pixel: F)
where
    F: FnMut(u32, u32),
{
    let mut edges: Vec<(delaunator::Point, delaunator::Point)> = Vec::new();
    let num_vertices = vertices.len();

    // Create the edge list
    for i in 0..num_vertices {
        let start = vertices[i].clone();
        let end = vertices[(i + 1) % num_vertices].clone();
        edges.push((start, end));
    }

    // Find the bounding box of the polygon
    let mut min_y = vertices[0].y;
    let mut max_y = vertices[0].y;
    for vertex in vertices {
        if vertex.y < min_y {
            min_y = vertex.y;
        }
        if vertex.y > max_y {
            max_y = vertex.y;
        }
    }

    // Scanline algorithm
    for y in min_y.floor() as i32..=max_y.ceil() as i32 {
        let mut intersections: Vec<f64> = Vec::new();

        // Find intersections of the current scanline with polygon edges
        for (start, end) in &edges {
            if (start.y <= y as f64 && end.y > y as f64)
                || (end.y <= y as f64 && start.y > y as f64)
            {
                let dx = end.x - start.x;
                let dy = end.y - start.y;
                let intersect_x = start.x + dx * ((y as f64 - start.y) / dy);
                intersections.push(intersect_x);
            }
        }

        // Sort intersections
        intersections.sort_by(|a, b| a.partial_cmp(b).unwrap());

        // Fill pixels between pairs of intersections
        for i in (0..intersections.len()).step_by(2) {
            if i + 1 < intersections.len() {
                let x_start = intersections[i];
                let x_end = intersections[i + 1];
                for x in x_start.ceil() as i32..=x_end.floor() as i32 {
                    if x >= 0 && y >= 0 {
                        draw_pixel(x as u32, y as u32);
                    }
                }
            }
        }
    }
}

fn centroid(vertices: &[delaunator::Point]) -> delaunator::Point {
    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut area = 0.0;

    let n = vertices.len();
    for i in 0..n {
        let x0 = vertices[i].x;
        let y0 = vertices[i].y;
        let x1 = vertices[(i + 1) % n].x;
        let y1 = vertices[(i + 1) % n].y;

        let a = x0 * y1 - x1 * y0;
        area += a;
        cx += (x0 + x1) * a;
        cy += (y0 + y1) * a;
    }

    area *= 0.5;
    cx /= (6.0 * area);
    cy /= (6.0 * area);

    delaunator::Point { x: cx, y: cy }
}

fn scale(vertices: &[delaunator::Point], scale_factor: f64) -> Vec<delaunator::Point> {
    let centroid = centroid(vertices);
    vertices
        .iter()
        .map(|v| delaunator::Point {
            x: centroid.x + scale_factor * (v.x - centroid.x),
            y: centroid.y + scale_factor * (v.y - centroid.y),
        })
        .collect()
}

fn scale_center(
    vertices: &[delaunator::Point],
    center: &delaunator::Point,
    scale_factor: f64,
) -> Vec<delaunator::Point> {
    vertices
        .iter()
        .map(|v| delaunator::Point {
            x: center.x + scale_factor * (v.x - center.x),
            y: center.y + scale_factor * (v.y - center.y),
        })
        .collect()
}

fn scale_floody(vertices: &[delaunator::Point], distance: f64) -> Vec<delaunator::Point> {
    let n = vertices.len();
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let x0 = Rector::from(&vertices[i]);
        let xn = Rector::from(&vertices[(i + 1) % n]);
        let xp = Rector::from(&vertices[if i == 0 { n - 1 } else { i - 1 }]);

        let vp = xp - x0;
        let vp = vp / vp.norm();

        let vn = xn - x0;
        let vn = vn / vn.norm();

        let nor = vp - vn;
        let nor = nor / nor.norm();

        let nor = Rector {
            x: -nor.y,
            y: nor.x,
        };

        let coss = vp * nor;
        let sinn = (1.0 - coss.powi(2)).sqrt();

        let (a, b, c) = (vp.norm(), vn.norm(), (xn - xp).norm());
        let incenter = (xp * a + x0 * b + xn * c) / (a + b + c);
        let d2in = (incenter - x0).norm();

        let s = d2in.min(distance / sinn);

        let outv = x0 + nor * s;

        out.push(delaunator::Point::from(outv));
    }
    out
}

fn bowyer_watson(x: u32, y: u32, points: Vec<Point2<f64>>) -> Vec<Rc<Triangle>> {
    let s = points.len() * 2 - 2;

    // let mut triangles: BTreeSet<Rc<Triangle>> = BTreeSet::new();
    let mut triangles: HashSet<Rc<Triangle>> = HashSet::with_capacity(s);
    // let mut triangles: Vec<Rc<Triangle>> = Vec::with_capacity(s);
    let mut bad_triangles = Vec::with_capacity(s);
    let mut polygon = Vec::with_capacity(s);

    let super_triangle = Triangle::new(
        Point2::new(-4.0 * x as f64, -4.0 * x as f64),
        Point2::new(4.0 * x as f64, -4.0 * x as f64),
        Point2::new(0.0, 4.0 * x as f64),
    );

    let super_triangle_vertices = super_triangle.vertices;

    // let mut triangles = vec![Rc::new(super_triangle)];
    triangles.insert(Rc::new(super_triangle));

    let input_pts: Vec<[f64; 2]> = points.iter().map(|p| [p.x, p.y]).collect();
    let (mut pts, bits) = hilbert::point_list::make_points_f64(&input_pts, 0, None, None, 10.0);
    hpoint::Point::hilbert_sort(&mut pts, bits);

    let sorted_pts: Vec<Point2<f64>> = pts
        .iter()
        .map(|p| {
            let [x, y] = p.get_coordinates()[..] else {
                unimplemented!();
            };
            Point2::new(x as f64 / 10.0, y as f64 / 10.0)
        })
        .collect();

    // let mut bad_triangles = Vec::new();
    for (point_index, point) in sorted_pts.iter().enumerate() {
        for triangle in triangles.iter() {
            if triangle.circumcircle_contains(&point) {
                get_region(triangle, &point, &mut bad_triangles);
                // eprintln!("Found bad triangle after {}", i);

                break;
            }
        }

        for triangle in &bad_triangles {
            for (edge_index, edge) in triangle.edges.iter().enumerate() {
                let is_shared = bad_triangles
                    .iter()
                    .filter(|t| t.edges.contains(edge))
                    .count()
                    > 1;
                if !is_shared {
                    polygon.push((
                        triangle.neighbours.borrow()[edge_index].upgrade(),
                        edge.clone(),
                    ));
                }
            }
        }

        for t in &bad_triangles {
            t.say_goodbye();
            triangles.remove(t);
        }

        for (edge_trian, edge) in &polygon {
            let trian = Rc::new(Triangle::new(edge.start, edge.end, *point));
            add_neighbours(&triangles, &trian);
            if let Some(et) = edge_trian {
                // eprintln!("triangle at poly edge");
                trian.update_neighbourhood(&et);
            }
            // triangles.push(trian);
            triangles.insert(trian);
        }

        bad_triangles.clear();
        polygon.clear();

        print_progress(point_index + 1, points.len());
    }

    triangles.retain(|t| {
        !super_triangle_vertices.contains(&t.vertices[0])
            && !super_triangle_vertices.contains(&t.vertices[1])
            && !super_triangle_vertices.contains(&t.vertices[2])
    });

    // triangles
    triangles.iter().cloned().collect()
}

fn generate_voronoi(x: u32, y: u32) -> image::ImageBuffer<image::Rgb<f32>, Vec<f32>> {
    let mut imgbuf = image::ImageBuffer::new(x, y);

    let mut rng = rand::thread_rng();
    // let points: Vec<Point2<f64>> = (0..1000)
    //     .map(|_| Point2::new(rng.gen_range(0.0..x as f64), rng.gen_range(0.0..y as f64)))
    //     .collect();

    let points: Vec<delaunator::Point> = (0..1000)
        .map(|_| delaunator::Point {
            x: rng.gen_range(0.0..x as f64),
            y: rng.gen_range(0.0..y as f64),
        })
        .collect();

    let start = Instant::now();
    // let triangles = bowyer_watson(x, y, points);
    let result = delaunator::triangulate(&points);

    eprint!("Bowyer watson took: {:?}", start.elapsed());
    for (_px, _py, pixel) in imgbuf.enumerate_pixels_mut() {
        *pixel = image::Rgb([0.5, 0.0, 0.0]);
    }

    for (i, triangle) in result.triangles.chunks(3).enumerate() {
        // eprintln!("--------------------------------------------------");
        // eprintln!("Triangle: {:?}", triangle);
        // for (j, n) in triangle.neighbours.borrow().iter().enumerate() {
        //     if let Some(nu) = n.upgrade() {
        //         eprintln!("===================");
        //         eprintln!("nachbar{:}: {:?}", j, nu);
        //     }
        // }

        let [i1, i2, i3] = triangle else {
            unimplemented!();
        };
        let t = [
            points[*i1].clone(),
            points[*i2].clone(),
            points[*i3].clone(),
        ];
        // let min_x = t.iter().map(|v| v.x).fold(f64::INFINITY, f64::min).floor() as u32;
        // let max_x = t
        //     .iter()
        //     .map(|v| v.x)
        //     .fold(f64::NEG_INFINITY, f64::max)
        //     .ceil() as u32;
        // let min_y = t.iter().map(|v| v.y).fold(f64::INFINITY, f64::min).floor() as u32;
        // let max_y = t
        //     .iter()
        //     .map(|v| v.y)
        //     .fold(f64::NEG_INFINITY, f64::max)
        //     .ceil() as u32;
        // let c = i as f32 / 10.0;
        //

        let temp_triangle = Triangle::new(
            Point2::new(t[0].x, t[0].y),
            Point2::new(t[1].x, t[1].y),
            Point2::new(t[2].x, t[2].y),
        );

        let c1 = temp_triangle.circumcenter.x as f32 / x as f32;
        let c2 = temp_triangle.circumcenter.y as f32 / y as f32;
        let c3 = smallest_angle((t[0].x, t[0].y), (t[1].x, t[1].y), (t[2].x, t[2].y));

        let edge_width = (1.0 - c1) * (1.0 - c2) * 10.0;
        let edge_width = edge_width.clamp(0.0, 10.0);
        let t2 = scale_floody(&t, edge_width as f64);

        scanline_fill(&t2, |a, b| {
            imgbuf.put_pixel(
                a.min(x),
                b.min(y),
                image::Rgb([c1, c2, 1.0 - c3 as f32 * 3.0 / 3.2]),
            )
        });

        // for x in min_x..=max_x {
        //     for y in min_y..=max_y {
        //         let point = Point2::new(x as f64, y as f64);
        //         if temp_triangle.contains_point(&point) {
        //             imgbuf.put_pixel(x, y, image::Rgb([c1, c2, 1.0 - c3 as f32 * 3.0 / 3.2]));
        //         }
        //     }
        // }
    }

    imgbuf
}

impl Noise for Voronoi {
    fn generate(&self, x: u32, y: u32) -> image::ImageBuffer<image::Rgb<f32>, Vec<f32>> {
        generate_voronoi(x, y)
    }
}

struct PolygonSet {
    vertices: Vec<usize>,
    caesures: Vec<usize>,
    halfedges: Vec<usize>,
}

fn to_voronoi(triangulation: delaunator::Triangulation) -> PolygonSet {
    let delaunator::Triangulation {
        triangles,
        halfedges,
        ..
    } = triangulation;

    let mut pts_to_triangle = Vec::with_capacity(triangles.len());
    for (i, p) in triangles.iter().enumerate() {
        if pts_to_triangle[*p] != 0 {
            pts_to_triangle[*p] = i;
        }
    }

    for (i, pt) in pts_to_triangle.iter().enumerate() {
        while true {
            temp_triangle = 
        }
    }
}

#[cfg(test)]
mod tests {

    use hilbert::point_list;

    use super::*;

    fn get_region2(start: &Rc<Triangle>) -> Vec<Weak<Triangle>> {
        let mut res = vec![Rc::downgrade(start)];

        fn inner(si: &Rc<Triangle>, r: &mut Vec<Weak<Triangle>>) {
            si.checked.replace(true);
            for neighbour_weak_ref in si.neighbours.borrow().iter() {
                if let Some(n_ref) = neighbour_weak_ref.upgrade() {
                    if !*n_ref.checked.borrow() {
                        r.push(Rc::downgrade(&n_ref));
                        inner(&n_ref, r);
                    }
                }
            }
        }

        inner(start, &mut res);
        res
    }

    fn print_trians<T: AsRef<[Rc<Triangle>]>>(triangles: T) {
        for t in triangles.as_ref() {
            println!("{}", t);
            println!("Neighbours");
            for n in t.neighbours.borrow().iter() {
                if let Some(nu) = n.upgrade() {
                    println!("{}", nu);
                }
            }
            println!("-----------------------");
        }
    }

    #[test]
    fn test_region() {
        let p1 = Point2::new(3.0, 0.0);
        let p2 = Point2::new(7.0, 0.0);
        let p3 = Point2::new(8.0, 3.0);
        let p4 = Point2::new(5.0, 5.0);
        let p5 = Point2::new(2.0, 3.0);
        let p6 = Point2::new(5.0, 2.0);

        let t1 = Rc::new(Triangle::new(p1, p2, p6));
        let t2 = Rc::new(Triangle::new(p2, p3, p6));
        let t3 = Rc::new(Triangle::new(p3, p4, p6));
        let t4 = Rc::new(Triangle::new(p4, p5, p6));
        let t5 = Rc::new(Triangle::new(p5, p1, p6));

        let triangles = vec![t1, t2, t3, t4, t5];
        add_neighbours(&triangles[1..].to_vec(), &triangles[0]);
        print_trians(&triangles);
        triangles[4].say_goodbye();

        println!("OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO");
        print_trians(&triangles);

        let r = get_region2(&triangles[0]);

        println!("OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO");
        for t in r {
            if let Some(tu) = t.upgrade() {
                println!("{}", tu);
            } else {
                self::panic!("AAAAAAAA");
            }
        }
    }
    #[test]
    fn test_hilbert() {
        let point_list = vec![
            Point2::new(3.0, 0.0),
            Point2::new(7.0, 0.0),
            Point2::new(8.0, 3.0),
            Point2::new(5.0, 5.0),
            Point2::new(2.0, 3.0),
            Point2::new(5.0, 2.0),
        ];

        let pts: Vec<[f64; 2]> = point_list.iter().map(|p| [p.x, p.y]).collect();
        let (mut points, bits) = hilbert::point_list::make_points_f64(&pts, 0, None, None, 10.0);
        print!("hilbert points 1==== {:?}", points);
        hpoint::Point::hilbert_sort(&mut points, bits);
        print!("hilbert points 2===={:?}", points);

        for p in &points {
            println!("hilbert point?????{:?}", p.hilbert_transform(bits));
            println!("hilbert point?????{:?}", p.get_coordinates());
        }
    }
}
