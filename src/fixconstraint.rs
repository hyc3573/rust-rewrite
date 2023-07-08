use crate::utils::*;
use crate::constraint::*;
use crate::point::*;

pub struct FixConstraint {
    pub cid: usize,
    pub pids: [usize; 1],
    pub x0: Real,
    pub y0: Real,
    pub z0: Real,
}

impl FixConstraint {
    pub fn init(&mut self, x0: Real, y0: Real, z0: Real) {
        self.x0 = x0;
        self.y0 = y0;
        self.z0 = z0;
    }
}

impl Constraint<1> for FixConstraint {
    fn new(pointid: [usize; 1], constraintid: usize) -> Self {
        FixConstraint {
            cid: constraintid,
            pids: pointid,
            x0: 0.,
            y0: 0.,
            z0: 0.,
        }
    }

    fn jacobian(
        &self,
        point: &[Point],
        j: &mut ns::coo::CooMatrix<Real>,
        jdot: &mut ns::coo::CooMatrix<Real>,
        c: &mut ns::coo::CooMatrix<Real>,
        cdot: &mut ns::coo::CooMatrix<Real>,
    ) {
        // C = 0.5((x-x0)^2 + (y-y0)^2 + (z-z0)^2)
        // J = (x-x0, y-y0, z-z0)

        let p = point[self.pids[0]];

        j.push(self.cid, 3 * p.id, p.x.x);
        j.push(self.cid, 3 * p.id + 1, p.x.y);
        j.push(self.cid, 3 * p.id + 2, p.x.z);

        // Jdot
        jdot.push(self.cid, 3 * p.id, p.v.x);
        jdot.push(self.cid, 3 * p.id + 1, p.v.y);
        jdot.push(self.cid, 3 * p.id + 1, p.v.z);

        c.push(self.cid, 0,
               0.5*(p.x.x-self.x0)*(p.x.x-self.x0) + 
               0.5*(p.x.y-self.y0)*(p.x.y-self.y0) + 
               0.5*(p.x.z-self.z0)*(p.x.z-self.z0)
        );

        cdot.push(self.cid, 0,
               p.x.x*p.v.x + 
               p.x.y*p.v.y + 
               p.x.z*p.v.z
        );
    }

    fn point_ids(&self) -> [usize; 1] {
        self.pids
    }
}
