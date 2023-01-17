#pragma once
#ifndef Scene_H
#define Scene_H

#include <vector>
#include <memory>
#include <string>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Rope;
class Particle;
class MatrixStack;
class Program;
class Shape;

class Scene
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	Scene();
	virtual ~Scene();
	
	void load(const std::string &RESOURCE_DIR);
	void init();
	void tare();
	void reset();
	void step();
	
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog) const;
	void drawSimple(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> progSimple) const;

	void setEndPos(Eigen::Vector3d endPosIn) { curRopeEndPos = endPosIn; }
	void setGrav(Eigen::Vector3d gravIn) { grav = gravIn; }
	void setRopeStiffness(double stiffIn);
	void tightenRopeSprings();
	void simSpheres(bool spheresOn) {showSpheres = spheresOn; }
	
	double getTime() const { return t; }
	
private:
	double t;
	double h;
	Eigen::Vector3d grav;

	Eigen::Vector3d curRopeEndPos;
	bool showSpheres;
	
	std::shared_ptr<Shape> sphereShape;
	std::shared_ptr<Rope> rope;
	std::vector< std::shared_ptr<Particle> > spheres;
};

#endif
