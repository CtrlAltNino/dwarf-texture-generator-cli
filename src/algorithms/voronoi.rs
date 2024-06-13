use clap::Args;
extern crate nalgebra as na;
extern crate rand;

use std::fmt;

use crate::algorithms::Noise;
use core::panic;
use na::{Point2, Vector2};
use rand::Rng;
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

impl PartialEq for Triangle {
    fn eq(&self, other: &Self) -> bool {
        self.vertices.iter().all(|v| other.vertices.contains(v))
    }
}

fn add_neighbours(triangles: &Vec<Rc<Triangle>>, new_triangle: &Rc<Triangle>) {
    for t in triangles {
        new_triangle.update_neighbourhood(t);
    }
}

fn get_region(start: &Rc<Triangle>, point: &Point2<f64>) -> Vec<Weak<Triangle>> {
    let mut res = vec![Rc::downgrade(start)];

    fn inner(si: &Rc<Triangle>, p: &Point2<f64>, r: &mut Vec<Weak<Triangle>>) {
        si.checked.replace(true);
        for neighbour_weak_ref in si.neighbours.borrow().iter() {
            if let Some(n_ref) = neighbour_weak_ref.upgrade() {
                if !*n_ref.checked.borrow() && n_ref.circumcircle_contains(p) {
                    r.push(Rc::downgrade(&n_ref));
                    inner(&n_ref, p, r);
                }
            }
        }
    }

    inner(start, point, &mut res);
    res
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

fn bowyer_watson(x: u32, y: u32, points: Vec<Point2<f64>>) -> Vec<Rc<Triangle>> {
    let super_triangle = Triangle::new(
        Point2::new(-4.0 * x as f64, -4.0 * x as f64),
        Point2::new(4.0 * x as f64, -4.0 * x as f64),
        Point2::new(0.0, 4.0 * x as f64),
    );

    let super_triangle_vertices = super_triangle.vertices;

    let mut triangles = vec![Rc::new(super_triangle)];
    let mut recently_added: Vec<Triangle> = vec![];

    for point in points {
        let mut bad_triangles_weak = Vec::new();
        for (i, triangle) in triangles.iter().enumerate() {
            if triangle.circumcircle_contains(&point) {
                bad_triangles_weak = get_region(triangle, &point);
                eprintln!("Found bad triangle after {}", i);

                break;
            }
        }

        for neighbour_weak_ref in &bad_triangles_weak {
            if let Some(n_ref) = neighbour_weak_ref.upgrade() {
                n_ref.checked.replace(false);
            } else {
                panic!("asdjfkasjdfas  aaaaaaaaaaaaaaaaaaaaaaaaa");
            }
        }
        // eprintln!("bad triangles{:}", bad_triangles_weak.len());

        let mut bad_triangles = Vec::new();

        for t_weak in bad_triangles_weak {
            if let Some(r) = t_weak.upgrade() {
                bad_triangles.push(r);
            } else {
                panic!("asdfasdfasdfasdf");
            }
        }

        let mut polygon = Vec::new();

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

        triangles.retain(|t| {
            if bad_triangles.contains(t) {
                t.say_goodbye();
                // eprintln!("saying good bye");
                false
            } else {
                // eprintln!("staying");
                true
            }
        });

        for (edge_trian, edge) in polygon {
            let trian = Rc::new(Triangle::new(edge.start, edge.end, point));
            add_neighbours(&triangles, &trian);
            if let Some(et) = edge_trian {
                // eprintln!("triangle at poly edge");
                trian.update_neighbourhood(&et);
            }
            triangles.push(trian);
        }
    }

    triangles.retain(|t| {
        !super_triangle_vertices.contains(&t.vertices[0])
            && !super_triangle_vertices.contains(&t.vertices[1])
            && !super_triangle_vertices.contains(&t.vertices[2])
    });

    triangles
}

fn generate_voronoi(x: u32, y: u32) -> image::ImageBuffer<image::Rgb<f32>, Vec<f32>> {
    let mut imgbuf = image::ImageBuffer::new(x, y);

    let mut rng = rand::thread_rng();
    let points: Vec<Point2<f64>> = (0..250)
        .map(|_| Point2::new(rng.gen_range(0.0..x as f64), rng.gen_range(0.0..y as f64)))
        .collect();

    let triangles = bowyer_watson(x, y, points);
    for (_px, _py, pixel) in imgbuf.enumerate_pixels_mut() {
        *pixel = image::Rgb([0.5, 0.0, 0.0]);
    }

    for (i, triangle) in triangles.iter().enumerate() {
        // eprintln!("--------------------------------------------------");
        // eprintln!("Triangle: {:?}", triangle);
        // for (j, n) in triangle.neighbours.borrow().iter().enumerate() {
        //     if let Some(nu) = n.upgrade() {
        //         eprintln!("===================");
        //         eprintln!("nachbar{:}: {:?}", j, nu);
        //     }
        // }
        let min_x = triangle
            .vertices
            .iter()
            .map(|v| v.x)
            .fold(f64::INFINITY, f64::min)
            .floor() as u32;
        let max_x = triangle
            .vertices
            .iter()
            .map(|v| v.x)
            .fold(f64::NEG_INFINITY, f64::max)
            .ceil() as u32;
        let min_y = triangle
            .vertices
            .iter()
            .map(|v| v.y)
            .fold(f64::INFINITY, f64::min)
            .floor() as u32;
        let max_y = triangle
            .vertices
            .iter()
            .map(|v| v.y)
            .fold(f64::NEG_INFINITY, f64::max)
            .ceil() as u32;
        let c = i as f32 / 10.0;

        let c1 = triangle.circumcenter.x as f32 / x as f32;
        let c2 = triangle.circumcenter.y as f32 / y as f32;

        for x in min_x..=max_x {
            for y in min_y..=max_y {
                let point = Point2::new(x as f64, y as f64);
                if triangle.contains_point(&point) {
                    imgbuf.put_pixel(x, y, image::Rgb([c1, c2, 1.0])); // Color the triangle red
                }
            }
        }
    }

    imgbuf
}

impl Noise for Voronoi {
    fn generate(&self, x: u32, y: u32) -> image::ImageBuffer<image::Rgb<f32>, Vec<f32>> {
        generate_voronoi(x, y)
    }
}

#[cfg(test)]
mod tests {

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
}
