extern crate nalgebra as na;
extern crate nalgebra_sparse as ns;
use core::cell::Cell;
use std::io;
use std::cell::RefCell;
use std::vec;

type Real = f32;

struct Idgen {
    cnt: Cell<usize>,
}

impl Idgen {
    fn new() -> Self {
        return Idgen { cnt: Cell::new(0) };
    }

    fn get(&self) -> usize {
        let temp = self.cnt.get();
        self.cnt.set(temp + 1);
        return temp;
    }

    fn size(&self) -> usize {
        return self.cnt.get();
    }
}

enum Mass {
    Mass(Real),
    Invmass(Real),
}

#[derive(Clone, Copy)]
struct Point {
    x: na::SVector<Real, 3>,
    v: na::SVector<Real, 3>,
    fsum: na::SVector<Real, 3>,
    w: Real,
    id: usize,
}

impl Point {
    fn new(mass: Mass, id: usize) -> Self {
        let w: Real;
        match mass {
            Mass::Mass(m) => w = 1. / m,
            Mass::Invmass(m) => w = m,
        }

        return Point {
            v: na::Vector3::new(0., 0., 0.),
            x: na::Vector3::new(0., 0., 0.),
            fsum: na::Vector3::new(0., 0., 0.),
            w,
            id,
        };
    }

    fn integrate(&mut self, dt: Real) -> () {
        self.v += self.fsum * self.w * dt;
        self.x += self.v * dt;
        self.fsum = na::zero();
    }

    fn add_force(&mut self, f: &na::Vector3<Real>) {
        self.fsum += f;
    }

    fn add_force_view<'a>(&mut self, f: &na::VectorView3<'a, Real, na::Dyn, na::Dyn>) {
        self.fsum += f;
    }

    fn set_force(&mut self, f: &na::Vector3<Real>) {
        self.fsum.copy_from(f);
    }
}

fn get_point(mass: Mass, ps: &mut vec::Vec::<Point>) {
        ps.push(Point::new(mass, ps.len()));
}

// Fixed -- 1
// On-the-line -- 1
// Linkage -- 2
// Angle -- 3
trait Constraint<const N: usize> {
    fn new(pointid: [usize; N], constraintid: usize) -> Self where Self: Sized;
    fn jacobian(
        &self,
        points: &vec::Vec::<Point>,
        j: &mut ns::coo::CooMatrix<Real>,
        jdot: &mut ns::coo::CooMatrix<Real>,
        c: &mut ns::coo::CooMatrix<Real>,
        cdot: &mut ns::coo::CooMatrix<Real>,
    );

    fn point_ids(&self) -> [usize; N];
    
    fn constraint_info(
        &self,
        point: &vec::Vec::<Point>,
        j: &mut ns::coo::CooMatrix<Real>,
        jdot: &mut ns::coo::CooMatrix<Real>,
        c: &mut ns::coo::CooMatrix<Real>,
        cdot: &mut ns::coo::CooMatrix<Real>,
        invmass: &mut ns::coo::CooMatrix<Real>,
        f: &mut ns::coo::CooMatrix<Real>,
        v: &mut ns::coo::CooMatrix<Real>,
    ) {
        for pid in self.point_ids() {
            let p = point[pid];
            invmass.push(3 * p.id, 3 * p.id, p.w);
            invmass.push(3 * p.id + 1, 3 * p.id + 1, p.w);
            invmass.push(3 * p.id + 2, 3 * p.id + 2, p.w);

            f.push(3 * p.id, 0, p.fsum.x);
            f.push(3 * p.id + 1, 0, p.fsum.y);
            f.push(3 * p.id + 2, 0, p.fsum.z);
            println!("{}", p.fsum);

            v.push(3 * p.id, 0, p.v.x);
            v.push(3 * p.id + 1, 0, p.v.y);
            v.push(3 * p.id + 2, 0, p.v.z);
        }

        self.jacobian(point, j, jdot, c, cdot);
    }
}

struct FixConstraint {
    cid: usize,
    pids: [usize; 1],
    x0: Real,
    y0: Real,
    z0: Real,
}

impl FixConstraint {
    fn init(&mut self, x0: Real, y0: Real, z0: Real) -> () {
        self.x0 = x0;
        self.y0 = y0;
        self.z0 = z0;
    }
}

impl Constraint<1> for FixConstraint {
    fn new(pointid: [usize; 1], constraintid: usize) -> Self {
        return FixConstraint {
            cid: constraintid,
            pids: pointid,
            x0: 0.,
            y0: 0.,
            z0: 0.,
        };
    }

    fn jacobian(
        &self,
        point: &vec::Vec::<Point>,
        j: &mut ns::coo::CooMatrix<Real>,
        jdot: &mut ns::coo::CooMatrix<Real>,
        c: &mut ns::coo::CooMatrix<Real>,
        cdot: &mut ns::coo::CooMatrix<Real>,
    ) {
        // C = 1/2((x-x0)^2 + (y-y0)^2 + (z-z0)^2)
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
               1/2*(p.x.x-x0)*(p.x.x-x0) + 
               1/2*(p.x.y-y0)*(p.x.y-y0) + 
               1/2*(p.x.z-z0)*(p.x.z-z0)
        );

        cdot.push(self.cid, 0,
               (p.x.x-x0)*(p.v.x) + 
               (p.x.y-y0)*(p.v.y) + 
               (p.x.z-z0)*(p.v.z)
        );
    }

    fn point_ids(&self) -> [usize; 1]{
        return self.pids;
    }
}

struct LinkConstraint {
    cid: usize,
    pids: [usize; 2],
    d: Real,
}

impl LinkConstraint{
    fn init(&mut self, d: Real) -> () {
        self.d = d;
    }
}

impl Constraint<2> for LinkConstraint {
    fn new(pointid: [usize; 2], constraintid: usize) -> Self {
        return LinkConstraint {
            cid: constraintid,
            pids: pointid,
            d: 0.,
        }
    }

    fn jacobian(
        &self,
        point: &vec::Vec::<Point>,
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
               1/2*(
                   (p1.x.x-p2.x.x)*(p1.x.x-p2.x.x) + 
                   (p1.x.y-p2.x.y)*(p1.x.y-p2.x.y) + 
                   (p1.x.z-p2.x.z)*(p1.x.z-p2.x.z) + 
                   - self.r*self.r
               )
        );

        cdot.push(self.cid, 0,
                  (p1.x.x-p2.x.x)*(p1.v.x-p2.v.x) + 
                  (p1.x.y-p2.x.y)*(p1.v.y-p2.v.y) + 
                  (p1.x.z-p2.x.z)*(p1.v.z-p2.v.z)
        );
    }

    fn point_ids(&self) -> [usize; 2] {
        return self.pids;
    }
}

struct ConstraintFactory<'a, const N: usize> {
    idgen: &'a Idgen,
}

impl<'a, const N: usize> ConstraintFactory<'a, N> {
    fn new(idgen: &Idgen) -> ConstraintFactory<N> {
        return ConstraintFactory {
            idgen,
        };
    }

    fn get<T: Constraint<N>>(&self, points: [Point; N]) -> T {
        return T::new(points.map(|p: Point| -> usize { p.id }), self.idgen.get());
    }

    fn size(&self) -> usize {
        return self.idgen.size();
    }
}

fn main() {
    println!("Hello, world!");

    let cidgen  = Idgen::new();
    let constraint1factory = ConstraintFactory::<1>::new(&cidgen);
    let constraint2factory = ConstraintFactory::<2>::new(&cidgen);

    // let mut p1 = pointfactory.get(Mass::Mass(1.));
    let mut ps = vec::Vec::<Point>::new();
    get_point(Mass::Mass(1.), &mut ps);
    ps[0].x += na::Vector3::new(2., 2., 2.);
    get_point(Mass::Mass(1.), &mut ps);
    ps[1].x += na::Vector3::new(7., 2., 2.);
    let mut c1 = constraint1factory.get::<FixConstraint>([ps[0]]);
    c1.init(0., 0., 0.);
    let mut c2 = constraint2factory.get::<LinkConstraint>([ps[0], ps[1]]);
    c2.init(5.);

    let c1s: vec::Vec::<Box<dyn Constraint<1>>> = vec!(
        Box::new(c1),
    );
    let c2s: vec::Vec::<Box<dyn Constraint<2>>> = vec!(
        Box::new(c2),
    );

    loop {
        // add external force
        ps[1].add_force(&na::Vector3::new(1., 1., 1.));

        // todo extract this to new function: from
        let mut j =
            ns::coo::CooMatrix::<Real>::new(constraint1factory.:)size(), ps.len() * 3);
        let mut jdot =
            ns::coo::CooMatrix::<Real>::new(constraint1factory.size(), ps.len() * 3);
        let mut invmass =
            ns::coo::CooMatrix::<Real>::new(ps.len() * 3, ps.len() * 3);
        let mut f = ns::coo::CooMatrix::<Real>::new(ps.len() * 3, 1);
        let mut v = ns::coo::CooMatrix::<Real>::new(ps.len() * 3, 1);
        let mut c = ns::coo::CooMatrix::<Real>::new(ps.len() * 3, 1);
        let mut cdot = ns::coo::CooMatrix::<Real>::new(ps.len() * 3, 1);

        for c in &c1s {
            c.constraint_info(&ps, &mut j, &mut jdot, &mut c, &mut cdot, &mut invmass, &mut f, &mut v);
        }
        for c in &c2s {
            c.constraint_info(&ps, &mut j, &mut jdot, &mut c, &mut cdot, &mut invmass, &mut f, &mut v);
        }

        let j = ns::convert::serial::convert_csr_dense(&ns::CsrMatrix::<Real>::from(&j));
        let jdot = ns::convert::serial::convert_csr_dense(&ns::CsrMatrix::<Real>::from(&jdot));
        let invmass =
            ns::convert::serial::convert_csr_dense(&ns::CsrMatrix::<Real>::from(&invmass));
        let f = ns::convert::serial::convert_csr_dense(&ns::CsrMatrix::<Real>::from(&f));
        let v = ns::convert::serial::convert_csr_dense(&ns::CsrMatrix::<Real>::from(&v));
        let c = ns::convert::serial::convert_csr_dense(&ns::CsrMatrix::<Real>::from(&v));
        let cdot = ns::convert::serial::convert_csr_dense(&ns::CsrMatrix::<Real>::from(&cdot));

        let left = &j * &invmass * &j.transpose();
        let right = -&jdot * &v - &j * &invmass * &f - 0.1*c - 0.1*cdot;
        println!("\njacobian: {}\njacobian_dot: {}\nleft: {}\nv: {}\nQ: {}\nright: {}", j, jdot, left, v, f, right);
        let lambda = left.try_inverse().unwrap()*right;
        let f = j.transpose()*lambda;
        println!("constraint force: {}", f);


        for i in 0..ps.len() {
            ps[i].add_force(&na::Vector3::new(f[3*i], f[3*i+1], f[3*i+2]));
            ps[i].integrate(0.1);
        }
        // to

        println!("resulting position: {}\nand\n{}", ps[0].x, ps[1].x);

        let mut asdf: String = String::new();
        asdf.reserve(2048);
        io::stdin().read_line(&mut asdf);
    }
}
