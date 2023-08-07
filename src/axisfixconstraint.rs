use crate::utils::*;
use crate::constraint::*;
use crate::point::*;

pub struct AxisFixConstraint<const N: usize> {
    pub cid: usize,
    pub pids: [usize; 1],
    pub x0: Real
}

impl<const N: usize> AxisFixConstraint<N> {
    pub fn init(&mut self, x0: Real) {
        self.x0 = x0;
    }
}

impl<const N: usize> Constraint<1> for AxisFixConstraint<{ N }> {
    fn new(pointid: [usize; 1], constraintid: usize) -> Self {
        Self {
            cid: constraintid,
            pids: pointid,
            x0: 0.0
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
        // C = x_N - x0
        // J = (1, 0, 0) or (0, 1, 0) or (0, 0, 1)

        let p = point[self.pids[0]];

        j.push(self.cid, 3*p.id+N, 1.0);

        c.push(self.cid, 0, p.x[N] - self.x0);
        cdot.push(self.cid, 0, p.v[N]);
    }

    fn point_ids(&self) -> [usize; 1] {
        self.pids
    }
}
