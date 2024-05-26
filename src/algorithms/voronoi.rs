use clap::Args;

use crate::algorithms::Noise;

use crate::composition::*;
use rand::prelude::*;

pub fn default_color(id: i32, center: Point, p: Point, size: Point) -> Rgb<f32> {
    Rgb([
        (center.x as f32 / size.x as f32),
        (center.y as f32 / size.y as f32),
        0.0,
    ])
}

pub fn default_color2(id: i32, center: Point, p: Point, size: Point) -> Rgb<f32> {
    let c = Rgb([
        (center.x as f32 / size.x as f32),
        (center.y as f32 / size.y as f32),
        0.0,
    ]);
    let d = center.distance(&p);
    c.interpolate(&Rgb([1.0, 1.0, 1.0]), d)
}

pub fn get_default_colorize() -> fn(i32, Point, Point, Point) -> Rgb<f32> {
    default_color2
}

#[derive(Args, Debug)]
pub struct Voronoi {
    // #[arg(short, long)]
    // seed: u32,
    #[arg(short, long)]
    cells: u32,

    #[arg(skip = get_default_colorize())]
    fgh: fn(i32, Point, Point, Point) -> Rgb<f32>,
}

extern crate image;

use image::{Rgb, Rgb32FImage};
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap};

// Define a point in 2D space
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct Point {
    x: i32,
    y: i32,
}

impl Point {
    fn distance(&self, other: &Point) -> f32 {
        let dx = (self.x - other.x) as f32;
        let dy = (self.y - other.y) as f32;
        (dx * dx + dy * dy).sqrt()
    }
}

// Event in the sweep line algorithm
#[derive(Clone, Copy)]
enum Event {
    Site(Point),
    Circle(Point, f32),
}

impl PartialEq for Event {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Event::Site(p1), Event::Site(p2)) => p1.eq(p2),
            (Event::Circle(p1, _), Event::Circle(p2, _)) => p1.eq(p2),
            _ => false,
        }
    }
}

impl Eq for Event {}

impl PartialOrd for Event {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Event {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self, other) {
            (Event::Site(p1), Event::Site(p2)) => p1.y.partial_cmp(&p2.y).unwrap(),
            (Event::Circle(_, y1), Event::Circle(_, y2)) => y1.partial_cmp(&y2).unwrap(),
            (Event::Site(_), Event::Circle(_, _)) => Ordering::Greater,
            (Event::Circle(_, _), Event::Site(_)) => Ordering::Less,
        }
    }
}

// Half-edge data structure for the Voronoi diagram
#[derive(Debug)]
struct HalfEdge {
    start: Point,
    end: Option<Point>,
    neighbor: Option<usize>,
}

// Event queue for the sweep line algorithm
struct EventQueue {
    events: BinaryHeap<Event>,
}

impl EventQueue {
    fn new() -> Self {
        EventQueue {
            events: BinaryHeap::new(),
        }
    }

    fn push(&mut self, event: Event) {
        self.events.push(event);
    }

    fn pop(&mut self) -> Option<Event> {
        self.events.pop()
    }
}

// Function to generate a Voronoi diagram
//
impl Voronoi {
    pub fn generate_voronoi(&self, x: u32, y: u32) -> Rgb32FImage {
        let mut image = Rgb32FImage::new(x, y);

        // Randomly generate sites
        let mut rng = rand::thread_rng();
        let mut sites = Vec::new();
        for _ in 0..self.cells {
            let site = Point {
                x: rng.gen_range(0..x as i32),
                y: rng.gen_range(0..y as i32),
            };
            sites.push(site);
        }

        // Event queue for the sweep line algorithm
        let mut event_queue = EventQueue::new();

        for site in &sites {
            event_queue.push(Event::Site(*site));
        }

        // Map from site to half-edge index
        let mut site_edge_map: HashMap<Point, usize> = HashMap::new();
        let mut edges = Vec::new();

        while let Some(event) = event_queue.pop() {
            match event {
                Event::Site(site) => {
                    // Handle site event
                    let edge_index = edges.len();
                    site_edge_map.insert(site, edge_index);
                    edges.push(HalfEdge {
                        start: site,
                        end: None,
                        neighbor: None,
                    });
                }
                Event::Circle(center, _radius) => {
                    // Handle circle event
                    if let Some(edge_index) = site_edge_map.get(&center) {
                        let mut edge = &mut edges[*edge_index];
                        edge.end = Some(center);
                    }
                }
            }
        }

        // Render Voronoi diagram
        for (px, py, pix) in image.enumerate_pixels_mut() {
            let query_point = Point {
                x: px as i32,
                y: py as i32,
            };
            let mut closest_site_index = 0;
            let mut min_distance = std::f32::MAX;
            for (i, site) in sites.iter().enumerate() {
                let distance = query_point.distance(site);
                if distance < min_distance {
                    min_distance = distance;
                    closest_site_index = i;
                }
            }
            let closest_site = &sites[closest_site_index];
            let color = (self.fgh)(
                closest_site_index as i32,
                *closest_site,
                query_point,
                Point {
                    x: x as i32,
                    y: y as i32,
                },
            );

            *pix = color;
        }

        image
    }
}

impl Noise for Voronoi {
    fn generate(&self, x: u32, y: u32) -> image::ImageBuffer<image::Rgb<f32>, Vec<f32>> {
        self.generate_voronoi(x, y)
    }
}
