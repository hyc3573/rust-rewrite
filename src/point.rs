use std::vec;
use nalgebra as na;

use crate::utils::*;

pub enum Mass {
    Mass(Real),
    Invmass(Real),
}

#[derive(Clone, Copy)]
pub struct Point {
    pub x: na::SVector<Real, 3>,
    pub v: na::SVector<Real, 3>,
    pub fsum: na::SVector<Real, 3>,
    pub w: Real,
    pub id: usize,
}

impl Point {
    pub fn new(mass: Mass, id: usize) -> Self {
        let w: Real = match mass {
            Mass::Mass(m) => 1. / m,
            Mass::Invmass(m) => m,
        };

        Point {
            v: na::Vector3::new(0., 0., 0.),
            x: na::Vector3::new(0., 0., 0.),
            fsum: na::Vector3::new(0., 0., 0.),
            w,
            id,
        }
    }

    pub fn integrate(&mut self, dt: Real) {
        self.v += self.fsum * self.w * dt;
        self.x += self.v * dt;
        self.fsum = na::zero();
    }

    pub fn add_force(&mut self, f: &na::Vector3<Real>) {
        self.fsum += f;
    }

    pub fn set_force(&mut self, f: &na::Vector3<Real>) {
        self.fsum.copy_from(f);
    }
}

pub fn get_point(mass: Mass, ps: &mut vec::Vec::<Point>) -> usize {
    ps.push(Point::new(mass, ps.len()));
    ps.len()-1
}
