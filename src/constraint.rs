use crate::point::*;
use crate::utils::*;
use crate::pv;

pub struct Solver {
    pub j: ns::coo::CooMatrix<Real>,
    pub jdot: ns::coo::CooMatrix<Real>,
    pub c: ns::coo::CooMatrix<Real>,
    pub cdot: ns::coo::CooMatrix<Real>,
    pub invmass: ns::coo::CooMatrix<Real>,
    pub f: ns::coo::CooMatrix<Real>,
    pub v: ns::coo::CooMatrix<Real>,
}

impl Solver {
    
    pub fn new(point_n: usize, constraint_n: usize) -> Self
    {
        type Coo = ns::coo::CooMatrix<Real>;
        Solver {
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
        let leftsize = left.shape().0;
        // let right = -&jdot * &v - &j * &invmass * &f - 0.0*&c - 0.0*&cdot;
        let right = -(&j*&invmass*&f);
        pv!(j, jdot, left, v, f, right, c, invmass);
        let lambda = left.try_inverse().unwrap_or(
            na::DMatrix::<Real>::zeros(leftsize, leftsize)
        )*right;
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
        solver: &mut Solver
    ) {
        let Solver { j, jdot, c, cdot, invmass, f, v } = solver;
        
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
