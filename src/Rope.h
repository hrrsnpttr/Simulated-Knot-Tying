#pragma once
#ifndef Rope_H
#define Rope_H

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>

class Particle;
class Spring;
class MatrixStack;
class Program;
class Shape;

class Rope
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	Rope(int verts,
		  const Eigen::Vector3d &xStart,
		  const Eigen::Vector3d &xEnd,
		  double mass,
		  double stiffness,
		  const std::shared_ptr<Shape> s);
	virtual ~Rope();
	
	void tare();
	void reset();
	void step(double h, const Eigen::Vector3d &grav, const std::vector< std::shared_ptr<Particle> > spheres, Eigen::Vector3d endPos);
	
	void drawParticles(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
	void drawConnectors(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> pSimple) const;

	void setStiffness(double stiffIn);
	void recreateSprings();
	
private:
	int verts;
	int n;
	std::vector< std::shared_ptr<Particle> > particles;
	std::vector< std::shared_ptr<Spring> > springs;
	
	Eigen::VectorXd v;
	Eigen::VectorXd f;
	Eigen::MatrixXd M;
	Eigen::MatrixXd K;
	
};

#endif
