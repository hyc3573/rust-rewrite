use nalgebra as na;
use std::vec;

use crate::{constraint::Solver, utils::*};

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

    pub fn point_info(&self, solver: &mut Solver) {
        let Solver { j, jdot, c, cdot, invmass, f, v } = solver;
        
        invmass.push(3 * self.id, 3 * self.id, self.w);
        invmass.push(3 * self.id + 1, 3 * self.id + 1, self.w);
        invmass.push(3 * self.id + 2, 3 * self.id + 2, self.w);

        f.push(3 * self.id, 0, self.fsum.x);
        f.push(3 * self.id + 1, 0, self.fsum.y);
        f.push(3 * self.id + 2, 0, self.fsum.z);

        v.push(3 * self.id, 0, self.v.x);
        v.push(3 * self.id + 1, 0, self.v.y);
        v.push(3 * self.id + 2, 0, self.v.z);
    }
}

pub fn get_point(mass: Mass, ps: &mut vec::Vec<Point>) {
    ps.push(Point::new(mass, ps.len()));
}
