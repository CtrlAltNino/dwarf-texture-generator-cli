use std::ops::{Add, Div, Mul, Sub};

use hilbert::normalize;

#[derive(Debug, Clone, Copy)]
pub struct Rector {
    pub x: f64,
    pub y: f64,
}

impl Rector {
    pub fn dot(&self, rhs: &Self) -> f64 {
        self.x * rhs.x + self.y * rhs.y
    }

    pub fn norm(&self) -> f64 {
        self.dot(self).sqrt()
    }

    pub fn normalize(&self) -> Self {
        *self / self.norm()
    }
}

trait Vector_Like {
    fn x(&self) -> f64;
    fn y(&self) -> f64;
}

impl Vector_Like for Rector {
    fn x(&self) -> f64 {
        let Rector { x, y } = self;
        *x
    }

    fn y(&self) -> f64 {
        let Rector { x, y } = self;
        *y
    }
}

impl From<&delaunator::Point> for Rector {
    fn from(value: &delaunator::Point) -> Self {
        let delaunator::Point { x, y } = value;
        Rector { x: *x, y: *y }
    }
}

impl From<Rector> for delaunator::Point {
    fn from(value: Rector) -> Self {
        let Rector { x, y } = value;
        delaunator::Point { x, y }
    }
}

impl Sub for Rector {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        Self::Output {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

impl Add for Rector {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::Output {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl Mul<f64> for Rector {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self::Output {
            x: self.x * rhs,
            y: self.y * rhs,
        }
    }
}

impl Div<f64> for Rector {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Self::Output {
            x: self.x / rhs,
            y: self.y / rhs,
        }
    }
}

impl Mul<Rector> for Rector {
    type Output = f64;

    fn mul(self, rhs: Self) -> f64 {
        self.dot(&rhs)
    }
}
