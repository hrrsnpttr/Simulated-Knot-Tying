#include <iostream>
#include <limits>

#include "Scene.h"
#include "Particle.h"
#include "Rope.h"
#include "Shape.h"
#include "Program.h"

using namespace std;
using namespace Eigen;

Scene::Scene() :
	t(0.0),
	h(1e-2),
	grav(0.0, 0.0, 0.0),
	showSpheres(false)
{
}

Scene::~Scene()
{
}

void Scene::load(const string &RESOURCE_DIR)
{
	
	sphereShape = make_shared<Shape>();
	sphereShape->loadMesh(RESOURCE_DIR + "sphere2.obj");
	
	// Units: meters, kilograms, seconds
	h = 5e-3; // 5e-3
	
	//grav << 0.0, -9.8, 0.0;
	//grav << 0.0, 0.0, 0.0;
	
	int verts = 70; // 70
	double mass = 0.1; // 0.1
	double stiffness = 1e2; // 1e1s 1e3g
	// Vector3d xStart(-0.25, 0.5, 0.0);
	// Vector3d xEnd(0.25, 0.5, 0.0);
	Vector3d xStart(-0.5, 0.5, 0.0);
	Vector3d xEnd(0.5, 0.5, 0.0);
	curRopeEndPos = xEnd;
	rope = make_shared<Rope>(verts, xStart, xEnd, mass, stiffness, sphereShape);
	
	auto sphere = make_shared<Particle>(sphereShape);
	spheres.push_back(sphere);
	sphere->r = 0.1;
	sphere->x = Vector3d(0.0, 0.475, 0.15);

}

void Scene::init()
{
	sphereShape->init();
}

void Scene::tare()
{
	for(int i = 0; i < (int)spheres.size(); ++i) {
		spheres[i]->tare();
	}
	rope->tare();
}

void Scene::reset()
{
	t = 0.0;
	for(int i = 0; i < (int)spheres.size(); ++i) {
		spheres[i]->reset();
	}
	Vector3d xEnd(0.5, 0.5, 0.0);
	curRopeEndPos = xEnd;
	rope->reset();
}

void Scene::setRopeStiffness(double stiffIn)
{
	rope->setStiffness(stiffIn);
}

void Scene::tightenRopeSprings() 
{
	rope->recreateSprings();
}

void Scene::step()
{
	t += h;
	
	// Move the sphere
	if(!spheres.empty()) {
		auto s = spheres.front();
		s->x(2) = 0.5 * sin(0.5*t);
	}
	
	// Simulate the rope
	if(showSpheres) {
		rope->step(h, grav, spheres, curRopeEndPos);
	} else {
		vector< shared_ptr<Particle> > spheresNull;
		rope->step(h, grav, spheresNull, curRopeEndPos);
	}

}

void Scene::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog) const
{
	glUniform3fv(prog->getUniform("kdFront"), 1, Vector3f(1.0, 1.0, 1.0).data());
	if(showSpheres) {
		for(int i = 0; i < (int)spheres.size(); ++i) {
			spheres[i]->draw(MV, prog);
		
		}
	}
	rope->drawParticles(MV, prog);
}

void Scene::drawSimple(shared_ptr<MatrixStack> MV, const shared_ptr<Program> progSimple) const
{
	//rope->drawConnectors(MV, progSimple);
}
