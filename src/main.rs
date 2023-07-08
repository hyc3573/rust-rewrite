extern crate nalgebra as na;
extern crate nalgebra_sparse as ns;
use core::cell::Cell;
use std::io;
use std::vec;

mod constraint;
use crate::constraint::*;

mod point;
use crate::point::*;

mod utils;
use crate::utils::*;

mod fixconstraint;
use crate::fixconstraint::*;

mod linkconstraint;
use crate::linkconstraint::*;

fn main() {
    let cidgen  = Idgen::new();
    let constraint1factory = ConstraintFactory::<1>::new(&cidgen);
    let constraint2factory = ConstraintFactory::<2>::new(&cidgen);

    // let mut p1 = pointfactory.get(Mass::Mass(1.));
    let mut ps = vec::Vec::<Point>::new();
    get_point(Mass::Mass(1.), &mut ps);
    ps[0].x += na::Vector3::new(2., 2., 2.);
    // get_point(Mass::Mass(1.), &mut ps);
    // ps[1].x += na::Vector3::new(5., 6., 2.);
    let mut c1 = constraint1factory.get::<FixConstraint>([ps[0]]);
    c1.init(2., 2., 2.);
    // let mut c2 = constraint2factory.get::<LinkConstraint>([ps[0], ps[1]]);
    // c2.init(5.);

    let c1s: vec::Vec::<Box<dyn Constraint<1>>> = vec!(
        Box::new(c1),
    );
    // let c2s: vec::Vec::<Box<dyn Constraint<2>>> = vec!(
    //     Box::new(c2),
    // );

    loop {
        // add external force
        ps[0].add_force(&na::Vector3::new(1., 0., 0.));
        // ps[1].add_force(&na::Vector3::new(1., 0., 0.));

        let mut solver = Solver::new(ps.len(), cidgen.size());

        for constraint in &c1s {
            constraint.constraint_info(&ps, &mut solver);
        }
        // for constraint in &c2s {
        //     constraint.constraint_info(&ps, &mut cinfo);
        // }

        for point in &ps {
            point.point_info(&mut solver);
        }

        let f = solver.calculate();
        pv!(f);


        for i in 0..ps.len() {
            ps[i].add_force(&na::Vector3::new(f[3*i], f[3*i+1], f[3*i+2]));
            ps[i].integrate(0.1);
        }

        println!("resulting position: {}", ps[0].x);

        let mut asdf: String = String::new();
        asdf.reserve(2048);
        io::stdin().read_line(&mut asdf).unwrap();
    }
}
