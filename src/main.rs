extern crate nalgebra as na;
extern crate nalgebra_sparse as ns;
use ::core::cell::Cell;
use std::io;
use std::vec;
use raylib::prelude::*;

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

mod axisfixconstraint;
use crate::axisfixconstraint::*;

fn main() {
    const WIDTH: Real = 640.0; const HEIGHT: Real = 480.0;
    let sin45: Real = (PI as Real/4.).sin();
    
    let (mut rl, thread) = raylib::init()
        .size(WIDTH as i32, HEIGHT as i32)
        .title("rust-rewrite")
        .build();
    
    let cidgen  = Idgen::new();
    let constraint1factory = ConstraintFactory::<1>::new(&cidgen);
    let constraint2factory = ConstraintFactory::<2>::new(&cidgen);

    // let mut p1 = pointfactory.get(Mass::Mass(1.));
    let mut ps = vec::Vec::<Point>::new();
    get_point(Mass::Mass(1.), &mut ps);
    ps[0].x += na::Vector3::new(WIDTH/2., HEIGHT/2., 0.);
    get_point(Mass::Mass(1.), &mut ps);
    ps[1].x += na::Vector3::new(WIDTH/2. + WIDTH/20.*sin45, HEIGHT/2. + WIDTH/20.*sin45, 0.);

    let mut c1 = constraint1factory.get::<AxisFixConstraint<0>>([ps[0]]);
    c1.init(WIDTH/2.);
    let mut c2 = constraint1factory.get::<AxisFixConstraint<1>>([ps[0]]);
    c2.init(HEIGHT/2.);
    let mut c3 = constraint1factory.get::<AxisFixConstraint<2>>([ps[0]]);
    c3.init(0.);

    let mut c4 = constraint2factory.get::<LinkConstraint>([ps[0], ps[1]]);
    c4.init(WIDTH/20.);

    let c1s: vec::Vec::<Box<dyn Constraint<1>>> = vec!(
        Box::new(c1),
        Box::new(c2),
        Box::new(c3),
    );
    let c2s: vec::Vec::<Box<dyn Constraint<2>>> = vec!(
        Box::new(c4),
    );

    let mut timer = std::time::Instant::now();
    while !rl.window_should_close() {
        let mut d = rl.begin_drawing(&thread);
        
        // add external force
        ps[0].add_force(&na::Vector3::new(0., WIDTH, 0.));
        ps[1].add_force(&na::Vector3::new(0., WIDTH, 0.));

        let mut solver = Solver::new(ps.len(), cidgen.size());

        for constraint in &c1s {
            constraint.constraint_info(&ps, &mut solver);
        }
        for constraint in &c2s {
            constraint.constraint_info(&ps, &mut solver);
        }

        for point in &ps {
            point.point_info(&mut solver);
        }

        let f = solver.calculate();
        // pv!(f);


        let dt = timer.elapsed().as_secs_f32();
        timer = std::time::Instant::now();
        for i in 0..ps.len() {
            ps[i].add_force(&na::Vector3::new(f[3*i], f[3*i+1], f[3*i+2]));
            ps[i].integrate(dt/10.0);
        }

        // println!("resulting position: {} and {}", ps[0].x, ps[1].x);
        // println!("distance {}", (ps[0].x - ps[1].x).norm());

        d.clear_background(Color::BLACK);

        for point in &ps {
            d.draw_circle(point.x.x as i32, point.x.y as i32, 5., Color::RED);
        }
        
        pv!((ps[0].x - ps[1].x).norm());
    }
}
