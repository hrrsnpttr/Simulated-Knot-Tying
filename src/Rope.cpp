#include <iostream>
#include <limits>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Rope.h"
#include "Particle.h"
#include "Spring.h"
#include "MatrixStack.h"
#include "Program.h"
#include "GLSL.h"

using namespace std;
using namespace Eigen;

shared_ptr<Spring> createSpring(const shared_ptr<Particle> p0, const shared_ptr<Particle> p1, double E) // Good
{
	auto s = make_shared<Spring>(p0, p1);
	s->E = E;
	Vector3d x0 = p0->x;
	Vector3d x1 = p1->x;
	Vector3d dx = x1 - x0;
	s->L = dx.norm();
	return s;
}

Rope::Rope(int verts,
		  const Eigen::Vector3d &xStart,
		  const Eigen::Vector3d &xEnd,
		  double mass,
		  double stiffness,
		  const shared_ptr<Shape> s) // Good
{
	assert(verts > 1);
	assert(mass > 0.0);
	assert(stiffness > 0.0);
	
	this->verts = verts;
	

	// Create particles

	this->n = 0; // size of global vector (do not count fixed vertices)
	double r = 0.01; // Used for collisions - 0.01
	int nVerts = verts;

	vector<Vector3d> pointsToFix{xStart, xEnd}; // Points to be fixed
	double massParticle = mass / nVerts; // Mass of each particle

	double xInc = (xEnd(0) - xStart(0)) / (nVerts - 1);
	double yInc = (xEnd(1) - xStart(1)) / (nVerts - 1);
	double zInc = (xEnd(2) - xStart(2)) / (nVerts - 1);

	for(int i = 0; i < nVerts; ++i) {
		auto p = make_shared<Particle>(s);
		particles.push_back(p);
		p->r = r;

		// Calculate world coords based on location relative to start and end points
		double newX = xStart(0) + (i * xInc);
		double newY = xStart(1) + (i * yInc);
		double newZ = xStart(2) + (i * zInc);

		// Store world coords in p->x
		Vector3d newCoords{newX, newY, newZ};
		p->x = newCoords;

		// Based on coords, determine if fixed
		bool fixPoint = false;
		for(int i = 0; i < pointsToFix.size(); i++) {
			if(pointsToFix[i] == newCoords) { fixPoint = true; break; }
		}

		// If fixed - Set 'fixed' to true, set 'i' to -1
		// If free  - Set 'fixed' to false, set 'i' to current 'n', increment 'n' by 3
		p->fixed = fixPoint;
		if(fixPoint) {
			p->i = -1;
		} else {
			p->i = n;
			n += 3;
		}

		// Store mass of particle in p->m
		p->m = massParticle;
	}


	// Create springs

	for(int i = 0; i < particles.size() - 1; ++i) {
		springs.push_back(createSpring(particles[i], particles[i + 1], stiffness)); 
	}


	// Allocate system matrices and vectors
	M.resize(n,n);
	K.resize(n,n);
	v.resize(n);
	f.resize(n);
}

Rope::~Rope() // Good
{
}

void Rope::recreateSprings() 
{
	double stiffness = springs[0]->E;
	springs.clear();
	for(int i = 0; i < particles.size() - 1; ++i) {
		springs.push_back(createSpring(particles[i], particles[i + 1], stiffness)); 
	}
}

void Rope::tare() // Good
{
	for(int k = 0; k < (int)particles.size(); ++k) {
		particles[k]->tare();
	}
}

void Rope::reset() // Good
{
	for(int k = 0; k < (int)particles.size(); ++k) {
		particles[k]->reset();
	}
}

void Rope::setStiffness(double stiffIn) 
{
	for(int i = 0; i < springs.size(); i++) {
		springs[i]->E = stiffIn;
	}
}

void Rope::step(double h, const Vector3d &grav, const vector< shared_ptr<Particle> > spheres, Vector3d endPos)
{
	M.setZero();
	K.setZero();
	v.setZero();
	f.setZero();

	Matrix3d I;
	I << 1, 0, 0,
	     0, 1, 0,
		 0, 0, 1;

	// Initial vector/matrix population
	for(int i = 0; i < particles.size(); i++) {

		int nLoc = particles[i]->i;

		// If not fixed
		if(nLoc != -1) {

			// Populate M matrix
			Matrix3d massTemp;
			double m = particles[i]->m;
			massTemp << m, 0, 0,
			            0, m, 0,
					    0, 0, m; 
			M.block<3,3>(nLoc, nLoc) = massTemp;

			// Populate v vector
			v.segment<3>(nLoc) = particles[i]->v;

			// Populate f vector
			f.segment<3>(nLoc) = (grav * m);
		}
	}

	// Apply spring forces
	for(int i = 0; i < springs.size(); i++) {
		
		// Calculate intermediate vars
		Vector3d dx = (springs[i]->p1)->x - (springs[i]->p0)->x; 
		double l = dx.norm();

		// Calculate Ks
		Matrix3d Ks = (springs[i]->E / (l * l)) * ((1 - ((l - springs[i]->L) / l)) * (dx * dx.transpose()) + ((l - springs[i]->L) / l) * (dx.dot(dx)) * I);  

		// Calculate fSpring
		Vector3d fSpring = springs[i]->E * (l - springs[i]->L) * (dx / l);

		// Apply fSpring to f vector
		if(springs[i]->p0->i != -1) {
			f.segment<3>(springs[i]->p0->i) += fSpring;
		}

		if(springs[i]->p1->i != -1) {
			f.segment<3>(springs[i]->p1->i) -= fSpring;
		}

		// Add Ks matrices to K
		if(springs[i]->p0->i != -1) {
			K.block<3, 3>(springs[i]->p0->i, springs[i]->p0->i) -= Ks;
		}

		if(springs[i]->p1->i != -1) {
			K.block<3, 3>(springs[i]->p1->i, springs[i]->p1->i) -= Ks;
		}

		if(springs[i]->p0->i != -1 && springs[i]->p1->i != -1) {
			// Off Diagonals
			K.block<3, 3>(springs[i]->p0->i, springs[i]->p1->i) += Ks;
			K.block<3, 3>(springs[i]->p1->i, springs[i]->p0->i) += Ks;
		}
	}

	// Handle Collisions
	double c = 100;
	double cFloor = 5;
	double cSelf = 50; // 30
	for(int i = 0; i < particles.size() - 1; i++) {
		
		// Sphere Collisions
		
		for(int j = 0; j < spheres.size(); j++) {
			
			// Calculate intermediate vars
			Vector3d dx = particles[i]->x - spheres[j]->x; 
			double l = dx.norm();
			Vector3d n = dx / l;
			double d = particles[i]->r + spheres[j]->r - l;

			// If collision detected
			if(d > 0) {

				// Apply collision force
				Vector3d fc = c * d * n;
				if(particles[i]->i != -1) {
					f.segment<3>(particles[i]->i) += fc;
				}

				// Alter K
				Matrix3d Kc = c * d * I;
				if(particles[i]->i != -1) {
					K.block<3, 3>(particles[i]->i, particles[i]->i) += Kc;
				}

			}
		}


		// Floor Collisions
		
		double dFloor = particles[i]->r - particles[i]->x(1);
		Vector3d nFloor;
		nFloor << 0.0, 1.0, 0.0;

		if(dFloor > 0) {
			
			// Apply collision force
			Vector3d fc = cFloor * dFloor * nFloor;
			if(particles[i]->i != -1) {
				f.segment<3>(particles[i]->i) += fc;
			}

			// Alter K
			Matrix3d Kc = cFloor * dFloor * I;
			if(particles[i]->i != -1) {
				K.block<3, 3>(particles[i]->i, particles[i]->i) += Kc;
			}
		}


		// Self Collisions

		int pad = 2; // Ignore this number of particles to the left and right of current
		
		for(int j = 0; j < particles.size(); j++) {
			
			// Skip particles within padded region
			if(j >= (i - pad) && j <= (i + pad)) { continue; }
			
			// Calculate intermediate vars
			Vector3d dx = particles[i]->x - particles[j]->x; 
			double l = dx.norm();
			Vector3d n = dx / l;
			double d = particles[i]->r + particles[j]->r - l;

			// If collision detected
			if(d > 0) {

				// Apply collision force
				Vector3d fc = cSelf * d * n;
				if(particles[i]->i != -1) {
					f.segment<3>(particles[i]->i) += fc;
				}

				// Alter K
				Matrix3d Kc = cSelf * d * I;
				if(particles[i]->i != -1) {
					K.block<3, 3>(particles[i]->i, particles[i]->i) += Kc;
				}

			}

		}

	}

	// Create A and b
	MatrixXd A = M - ((h * h) * K);
	VectorXd b = (M * v) + (h * f);

	// Solve Ax = b
	VectorXd x = A.ldlt().solve(b);

	// Update particle positions and velocities
	for(int i = 0; i < particles.size(); i++) {
		if(particles[i]->i != -1) {
			particles[i]->v = x.segment<3>(particles[i]->i);
			particles[i]->x += particles[i]->v * h;
			particles[i]->v = particles[i]->v * 0.99; // Artifical energy loss, mimics friction
		}
		if(i == particles.size() - 1) {
			particles[i]->x = endPos;
		}
	}

}

void Rope::drawParticles(shared_ptr<MatrixStack> MV, const shared_ptr<Program> p) const // Good
{
	// Draw Particles
	glUniform3fv(p->getUniform("kdFront"), 1, Vector3f(1.0, 1.0, 1.0).data());
	for(int i = 0; i < (int)particles.size(); ++i) {
		particles[i]->draw(MV, p);
	}
}

void Rope::drawConnectors(shared_ptr<MatrixStack> MV, const shared_ptr<Program> pSimple) const // Good
{
	glColor3f(0.5f, 0.0f, 0.0f);
	glLineWidth(4.0f); // 4
	glBegin(GL_LINE_STRIP);

	for(int i = 0; i < (int)particles.size(); ++i) {
		glVertex3f(particles[i]->x(0), particles[i]->x(1), particles[i]->x(2));
	}

	glEnd();
}
