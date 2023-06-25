use crate::point::*;
use crate::utils::*;
use crate::pv;

pub struct ConstraintInfo {
    j: ns::coo::CooMatrix<Real>,
    jdot: ns::coo::CooMatrix<Real>,
    c: ns::coo::CooMatrix<Real>,
    cdot: ns::coo::CooMatrix<Real>,
    invmass: ns::coo::CooMatrix<Real>,
    f: ns::coo::CooMatrix<Real>,
    v: ns::coo::CooMatrix<Real>,
}

impl ConstraintInfo {
    
    pub fn new(point_n: usize, constraint_n: usize) -> Self
    {
        type Coo = ns::coo::CooMatrix<Real>;
        ConstraintInfo {
            j: Coo::new(constraint_n, point_n*3),
            jdot: Coo::new(constraint_n, point_n*3),
            c: Coo::new(constraint_n, 1),
            cdot: Coo::new(constraint_n, 1),
            invmass: Coo::new(point_n*3, point_n*3),
            f: Coo::new(point_n*3, 1),
            v: Coo::new(point_n*3, 1),
        }
    }

    pub fn calculate(&self) -> na::DMatrix<Real> {
        let j = ns::convert::serial::convert_csr_dense(&ns::CsrMatrix::<Real>::from(&self.j));
        let jdot =
            ns::convert::serial::convert_csr_dense(&ns::CsrMatrix::<Real>::from(&self.jdot));
        let invmass =
            ns::convert::serial::convert_csr_dense(&ns::CsrMatrix::<Real>::from(&self.invmass));
        let f = ns::convert::serial::convert_csr_dense(&ns::CsrMatrix::<Real>::from(&self.f));
        let v = ns::convert::serial::convert_csr_dense(&ns::CsrMatrix::<Real>::from(&self.v));
        let c = ns::convert::serial::convert_csr_dense(&ns::CsrMatrix::<Real>::from(&self.c));
        let cdot = ns::convert::serial::convert_csr_dense(&ns::CsrMatrix::<Real>::from(&self.cdot));

        let left = &j * &invmass * &j.transpose();
        let right = -&jdot * &v - &j * &invmass * &f - 0.1*&c - 0.1*&cdot;
        pv!(j, jdot, left, v, f, right, c);
        let lambda = left.try_inverse().unwrap()*right;
        j.transpose()*lambda
    }
}

pub trait Constraint<const N: usize> {
    fn new(pointid: [usize; N], constraintid: usize) -> Self where Self: Sized;
    fn jacobian(
        &self,
        points: &[Point],
        j: &mut ns::coo::CooMatrix<Real>,
        jdot: &mut ns::coo::CooMatrix<Real>,
        c: &mut ns::coo::CooMatrix<Real>,
        cdot: &mut ns::coo::CooMatrix<Real>,
    );

    fn point_ids(&self) -> [usize; N];
    
    fn constraint_info(
        &self,
        point: &[Point],
        cinfo: &mut ConstraintInfo
    ) {
        let ConstraintInfo { j, jdot, c, cdot, invmass, f, v } = cinfo;
        
        for pid in self.point_ids() {
            let p = point[pid];
            invmass.push(3 * p.id, 3 * p.id, p.w);
            invmass.push(3 * p.id + 1, 3 * p.id + 1, p.w);
            invmass.push(3 * p.id + 2, 3 * p.id + 2, p.w);

            f.push(3 * p.id, 0, p.fsum.x);
            f.push(3 * p.id + 1, 0, p.fsum.y);
            f.push(3 * p.id + 2, 0, p.fsum.z);
            pv!(p.fsum);

            v.push(3 * p.id, 0, p.v.x);
            v.push(3 * p.id + 1, 0, p.v.y);
            v.push(3 * p.id + 2, 0, p.v.z);
        }

        self.jacobian(point, j, jdot, c, cdot);
    }
}

pub struct ConstraintFactory<'a, const N: usize> {
    idgen: &'a Idgen,
}

impl<'a, const N: usize> ConstraintFactory<'a, N> {
    pub fn new(idgen: &Idgen) -> ConstraintFactory<N> {
        ConstraintFactory {
            idgen,
        }
    }

    pub fn get<T: Constraint<N>>(&self, points: [Point; N]) -> T {
        T::new(points.map(|p: Point| -> usize { p.id }), self.idgen.get())
    }
}
