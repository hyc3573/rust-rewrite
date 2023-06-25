use crate::utils::*;
use crate::constraint::*;
use crate::point::*;

pub struct LinkConstraint {
    pub cid: usize,
    pub pids: [usize; 2],
    pub d: Real,
}

impl LinkConstraint{
    pub fn init(&mut self, d: Real) {
        self.d = d;
    }
}

impl Constraint<2> for LinkConstraint {
    fn new(pointid: [usize; 2], constraintid: usize) -> Self {
        LinkConstraint {
            cid: constraintid,
            pids: pointid,
            d: 0.,
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
        let p1 = point[self.pids[0]];
        let p2 = point[self.pids[1]];

        j.push(self.cid, 3*p1.id, p1.x.x - p2.x.x);
        j.push(self.cid, 3*p1.id + 1, p1.x.y - p2.x.y);
        j.push(self.cid, 3*p1.id + 2, p1.x.z - p2.x.z);

        j.push(self.cid, 3*p2.id, p2.x.x - p1.x.x);
        j.push(self.cid, 3*p2.id + 1, p2.x.y - p1.x.y);
        j.push(self.cid, 3*p2.id + 2, p2.x.z - p1.x.z);

        jdot.push(self.cid, 3*p1.id, p1.v.x - p2.v.x);
        jdot.push(self.cid, 3*p1.id + 1, p1.v.y - p2.v.y);
        jdot.push(self.cid, 3*p1.id + 2, p1.v.z - p2.v.z);

        jdot.push(self.cid, 3*p2.id, p2.v.x - p1.v.x);
        jdot.push(self.cid, 3*p2.id + 1, p2.v.y - p1.v.y);
        jdot.push(self.cid, 3*p2.id + 2, p2.v.z - p1.v.z);

        c.push(self.cid, 0,
               0.5*(
                   (p1.x.x-p2.x.x)*(p1.x.x-p2.x.x) + 
                   (p1.x.y-p2.x.y)*(p1.x.y-p2.x.y) + 
                   (p1.x.z-p2.x.z)*(p1.x.z-p2.x.z) + 
                   - self.d*self.d
               )
        );

        cdot.push(self.cid, 0,
                  (p1.x.x-p2.x.x)*(p1.v.x-p2.v.x) + 
                  (p1.x.y-p2.x.y)*(p1.v.y-p2.v.y) + 
                  (p1.x.z-p2.x.z)*(p1.v.z-p2.v.z)
        );
    }

    fn point_ids(&self) -> [usize; 2] {
        self.pids
    }
}
