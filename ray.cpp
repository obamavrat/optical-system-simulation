#include "vector.cpp"

#include <string>
#include <thread>

const int Num =128;



class ray {
public:
	vector T;
	vector R;
	bool activ = true;


	void set(vector T, vector R, double n) {
		this->activ = true;
		this->R = R;
		this->T = (T * n / T.Norm());
	}
	ray operator=(ray A) {
		this->T = A.T;
		this->R = A.R;
		this->activ = A.activ;
	}
	vector translate(double t) {
		vector res = this->R + this->T * t;
		return res;
	}
	void printT() {
		std::cout << std::to_string(this->T.x) << '\t' << std::to_string(this->T.y) << '\t' << std::to_string(this->T.z) << std::endl;
	}
	void printR() {
		std::cout << std::to_string(this->R.x) << '\t' << std::to_string(this->R.y) << '\t' << std::to_string(this->R.z) << std::endl;
	}
};

class plane {
public:
	vector a;
	vector b;
	vector r;
	void set(vector a, vector b, vector r) {
		/*if (a * b == a.Norm() * b.Norm()) {
			//тут надо сделать ошибку
		}*/
		this->r = r;
		if (a.Norm() != 1)  a.normalize();
		if (b.Norm() != 1)  b.normalize();
		if (a * b != 0) {
			b = b - a * (a * b) / (a * a);
			this->a = a;
			this->b = b;
		}
		else {
			this->a = a;
			this->b = b;
		}
	}
	vector normal() {
		vector res;
		res = this->a.X(this->b);
		res.normalize();
		return res;
	}
	vector transvers(ray L, bool& error) {
		vector Y = L.R - this->r;
		error = false;
		vector n = this->normal();
		double det = n * L.T;
		if (det == 0) { 
			error = true;
			return L.R;
		};
		double t = (Y * n) / (L.T * n);
		vector R = L.R - L.T * t;
		return R;
	}
	ray translate(ray L) {
		bool err = false;
		vector newR = this->transvers(L,err);
		if (err) { return L; }
		ray res;
		res.set(L.T, newR, L.T.Norm());
		return res;
	}
};
class cone {
public:
	double alpha;
	double z;
	void set(double alpha, double z) {
		this->alpha = alpha;
		this->z = z;
	}
	vector transvers(ray L, bool& error, double& d) {
		error = false;
		vector res;
		double a = L.T.x * L.T.x + L.T.y * L.T.y - abs(this->alpha) * L.T.z * L.T.z;
		double b = L.T.x * L.R.x * 2 + L.T.y * L.R.y * 2 - 2 * abs(this->alpha) * L.T.z * (L.R.z - this->z);
		double c = L.R.x * L.R.x + L.R.y * L.R.y - abs(this->alpha) * (L.R.z - this->z) * (L.R.z - this->z);
		double D = b * b - 4 * a * c;
		//std::cout << std::to_string(D) << std::endl;
		//std::cout << std::to_string(this->alpha) << std::endl;
		if (D < 0) {
			error = true;
			res.zero();
			return res;
		}
		else {
			vector res1;
			vector res2;
			double d1 = (-b + sqrt(D)) / (2 * a);
			double d2 = (-b - sqrt(D)) / (2 * a);
			res1 = L.T * d1 + L.R;
			res2 = L.T * d2 + L.R;
			//std::cout << std::to_string(L.T.z) << std::endl;
			//std::cout << std::to_string(res2.z) << '\t' << std::to_string(res1.z) << std::endl;
			//std::cout << '\t' << std::to_string(d2) << '\t' << std::to_string(d1) << std::endl;
			if (this->alpha > 0) {
				/*std::cout << std::to_string(res1.z)<<'\t' << std::to_string(this->z) << std::endl;
				std::cout << std::to_string(res1.z) << '\t' << std::to_string(this->z) << std::endl;*/
				if (res1.z > this->z && res2.z > this->z) {
					if (d1 > d2) { 
						d = d2;
						return res2; 
					}
					else { 
						d = d1;
						return res1;
					}
				}
				if (res1.z < this->z && res2.z < this->z) {
					error = true;
					//std::cout << "chtoto jsmislennoe" << std::endl;
					res.zero();
					return res;
				}
				if (res1.z < this->z && res2.z > this->z) {
					d = d2;
					return res2;
				}
				if (res1.z > this->z && res2.z < this->z) {
					d = d1;
					return res1;
				}
			}
			if (this->alpha < 0) {
				//std::cout << "chtoto jsmislennoe" << std::endl;
				//std::cout << std::to_string(res1.z - this->z) << '\t' << std::to_string(res2.z - this->z) << std::endl;
				if (res1.z > this->z && res2.z > this->z) {
					//std::cout << std::to_string(res2.z) << std::endl;
					error = true;
					res.zero();
				}
				if (res1.z < this->z && res2.z < this->z) {
					if (d1 > d2) {
						d = d2;
						return res2;
					}
					else {
						d = d1;
						return res1;
					}
				}
				if (res1.z < this->z && res2.z > this->z) {
					d = d1;
					return res1;
				}
				if (res1.z > this->z && res2.z < this->z) {
					d = d2;
					return res2;
				}
			}
			if (this->alpha == 0) {
				d = d1;
				return res1;
			}
		}
	}
	vector normal(vector2 r) {
		vector res;
		res.set(r.x, r.y, 0);
		res.z = -sqrt(abs(this->alpha)) / this->alpha * r.Norm();
		res.normalize();
		return res;
	}

};
class AxisSurface {
public:
	double z;
	double R[Num];
	double t[Num];
	void setGrid(double minR, double maxR) {
		double delta = (maxR - minR) / (Num - 1);
		for (int i = 0;i < Num;i++) this->R[i] = minR + i * delta;	
	}
	void setZ(double z) { this->z = z; }
	void setT(double* t) {
		for (int i = 0;i < Num;i++) {
			this->t[i] = t[i];
		}
	}
	vector transvers(ray L, bool& error,int& numbersubSurface) {
		bool errorTransCone[Num - 1];
		bool transSurf[Num - 1];
		double range[Num - 1];
		double alpha[Num - 1];
		double z0[Num - 1];
		double dR[Num - 1];
		double dt[Num - 1];
		vector a[Num - 1];
		cone tempCone;
		for (int i = 0;i < Num - 1;i++) {
			transSurf[i] = false;
			dR[i] = this->R[i + 1] - this->R[i];
			dt[i] = this->t[i + 1] - this->t[i];
			alpha[i] = pow(dR[i] / dt[i], 3) / abs(dR[i] / dt[i]);
			
			z0[i] = -this->R[i] * dt[i] / dR[i] + this->z + this->t[i];

			//std::cout << std::to_string(alpha[i]) << std::endl;

			tempCone.set(alpha[i], z0[i]);
			a[i] = tempCone.transvers(L, errorTransCone[i], range[i]);
			//std::cout << std::to_string(alpha[i]) << std::endl;
			//std::cout << std::to_string(a[i].x) << '\t' << std::to_string(a[i].y) << '\t' << std::to_string(a[i].z) << std::endl;
			if (!errorTransCone[i]) {
				//std::cout << "chtoto jsmislennoe" << std::endl;
				if (this->t[i + 1] > this->t[i]) {
					//std::cout << "chtoto jsmislennoe" << std::endl;
					if (this->t[i] + this->z < a[i].z && a[i].z < this->t[i + 1] + this->z) {
						transSurf[i] = true;
					}
				}
				if (this->t[i + 1] < this->t[i]) {
					//std::cout << "chtoto jsmislennoe" << std::endl;
					if (this->t[i + 1] + this->z < a[i].z && a[i].z < this->t[i] + this->z) {
						//std::cout << "chtoto jsmislennoe" << std::endl;
						transSurf[i] = true;
					}
				}
			}
		}
		double distanceToTrans = 0;
		int numberOfPoint = 0;
		bool startAlg = false;
		for (int i = 0;i < Num - 1;i++) {
			if (transSurf[i]) {
				if (!startAlg) {
					distanceToTrans = range[i];
					startAlg = true;
					numberOfPoint = i;
					continue;
				}
				if (distanceToTrans > range[i]) {
					distanceToTrans = range[i];
					numberOfPoint = i;
				}
			}
		}
		//std::cout<<std::to_string(z0[1]) << std::endl;
		if (!startAlg) {
			L.activ = false;
			a[0].zero();
			//std::cout << "dsgdsfg" << std::endl;
			return a[0];
		}
		else {
			numbersubSurface = numberOfPoint;
			error = false;
			return a[numberOfPoint];
		}
	}
		
	vector normal(vector r,int i) {
		cone P;
		double dR = this->R[i] - this->R[i + 1];
		double dt = this->t[i] - this->t[i + 1];
		double alpha = pow(dR / dt, 3) / abs(dR / dt);
		double z0 = -this->R[i] * dt / dR + this->z + this->t[i];
		P.set(alpha, z0);
		vector res = P.normal(r.proectionXY());
		return res;
	}
	ray mirrorReflectivity(ray L) {
		bool e = true;
		int i = 0;
		vector point = this->transvers(L, e, i);
		//std::cout << std::to_string(point.x) << '\t' << std::to_string(point.y) << '\t' << std::to_string(point.z) << std::endl;
		if (e) {
			//std::cout << "dfgfds" << std::endl;
			ray newL;
			newL.set(L.T, L.R, L.T.Norm());
			newL.activ = false;
			return newL;
		}
		else {
			vector N = this->normal(point, i);
			double a = L.T * N;
			vector newT;
			newT = (L.T + N * a * 2) * (-1);
			ray res;
			res.set(newT, point, L.T.Norm());
			return res;
		}

	}
	
	
};
class axisSpherSurface {
public:
	double R;
	double d;
	vector shift;
	void set(double R, double d) {
		this->R = R;
		this->d =  d;
		//std::cout << std::to_string(this->d) << std::endl;
		this->shift.set(0, 0, R + d );
	}
	vector center() {
		vector res;
		res.set(0, 0, this->d + this->R);
		return res;
	}
	vector transvers(ray L,bool& err) {
		//L.printR();
		vector center = this->shift;
		double D = ((L.R - center) * L.T) * ((L.R - center) * L.T) - ((L.R - center) * (L.R - center) - this->R * this->R);
		//std::cout << std::to_string(D)  << std::endl;
		if (D < 0) {
			
			err = true;
			return L.R;
		}
		double t = (center - L.R) * L.T + -abs(this->R) / this->R * sqrt(D);
		vector res = L.translate(t);
		return res;
	}
	vector normal(vector a) {
		vector res = (a - this->shift);
		//std::cout << std::to_string((a + this->shift).x) << "\t" << std::to_string((a + this->shift).y) << "\t" << std::to_string((a + this->shift).z) << std::endl;
		return res;
	}
	ray mirrorReflectivity(ray L) {
		bool error = false;
		
		vector trans = this->transvers(L,error);
		
		if (error) {
			std::cout << "dfgfd" << std::endl;
			ray res;
			res = L;
			res.activ = false;
			return L;
		}
		vector normal = this->normal(trans);
		
		normal.normalize();
		//normal.print();
		//std::cout << std::to_string(normal.z) << "\t" << std::to_string(normal.x) << std::endl;
		vector newT = (L.T - normal * (L.T * normal) * 2);
		ray res;
		
		//std::cout << std::to_string(newT.z) << "\t" << std::to_string(newT.x) << std::endl;
		
		res.set(newT, trans, L.T.Norm());
		
		return res;
	}
	ray refractive(ray L, double n1, double n2, bool& PVO) {
		bool err = false;
		vector trans = this->transvers(L, err);
		if (err) {
			ray res;
			res = L;
			res.activ = false;
			return res;
		}
		
		vector normal = this->normal(trans);
		vector su;
		su.setUnitVector(0, acos(0) );
		vector sv;
		sv.setUnitVector(acos(0) , acos(0) );
		vector Su = normal.X(su);
		vector Sv = normal.X(sv);
		
		double b = ((L.T * Su) * (Su * Sv) - (L.T * Sv) * (Su * Su)) / ((Su * Sv) * (Su * Sv) - (Su * Su) * (Sv * Sv));
		double a = L.T * Su / (Su * Su) - b * (Su * Sv) / (Su * Su);
		//std::cout << std::to_string(a)<<'\t' << std::to_string(a) << std::endl;
		vector newK = Su * a + Sv * b;
		//std::cout << std::to_string(newK*normal) << std::endl;
		double alphaSqr = (n2 * n2 - newK.sqrNorm()) / normal.sqrNorm();
		//std::cout << std::to_string(alphaSqr) << std::endl;
		if (alphaSqr < 0) {
			PVO = true;
			ray res = this ->mirrorReflectivity(L);
			return res;
		}
		
		newK = newK + normal * sqrt(alphaSqr);

		//std::cout << std::to_string(newK.sqrNorm()) << std::endl;
		if (newK * L.T < 0) {
			newK = newK - normal * sqrt(alphaSqr) * 2;
			//std::cout <<"dfjg" << std::endl;
		}
		ray res;
		res.set(newK, trans, n2);
		//std::cout << std::to_string(newK.x) << "\t" << std::to_string(newK.y) << "\t" << std::to_string(newK.z) << std::endl;
		//std::cout << std::to_string(newK.x) << "\t" << std::to_string(newK.y) << std::endl;
		return res;
	}
};
class emiter {

public:
	static const int numberRay = 12;
	ray E[numberRay][numberRay];
	void set(emiter e) {
		for (int i = 0;i < numberRay;i++) {
			for (int j = 0;j < numberRay;j++) {
				this->E[i][j] = e.E[i][j];
			}
		}
	}
	emiter reflectToMirror(AxisSurface S) {
		emiter res;
		ray v1;
		ray v2;
		ray v3;
		ray v4;
		ray V1;
		ray V2;
		ray V3;
		ray V4;
		int reduceNumberRay = numberRay / 4;
		for (int i = 0; i < numberRay;i++) {
			for (int j = 0; j < reduceNumberRay;j++) {
				V1 = this->E[i][j * 4];
				V2 = this->E[i][j * 4 + 1];
				V3 = this->E[i][j * 4 + 2];
				V4 = this->E[i][j * 4 + 3];
				std::thread th1([&v1, &V1, &S]() {
					if (V1.activ) {
						v1 = S.mirrorReflectivity(V1);
						//std::cout << std::to_string(res.E[i][j].T.x) << '\t' << std::to_string(res.E[i][j].T.y) << '\t' << std::to_string(res.E[i][j].T.z) << std::endl;
					}
					else {
						v1 = V1;
						//std::cout << std::to_string(i * 10 + j) << std::endl;
					}});
				std::thread th2([&v2, &V2, &S]() {
					if (V2.activ) {
						v2 = S.mirrorReflectivity(V2);
							//std::cout << std::to_string(res.E[i][j].T.x) << '\t' << std::to_string(res.E[i][j].T.y) << '\t' << std::to_string(res.E[i][j].T.z) << std::endl;
					}
					else {
						v2 = V2;
							//std::cout << std::to_string(i * 10 + j) << std::endl;
					}});
				std::thread th3([&v3, &V3, &S]() {
					if (V3.activ) {
						v3 = S.mirrorReflectivity(V3);
								//std::cout << std::to_string(res.E[i][j].T.x) << '\t' << std::to_string(res.E[i][j].T.y) << '\t' << std::to_string(res.E[i][j].T.z) << std::endl;
					}
					else {
						v3 = V3;
							//std::cout << std::to_string(i * 10 + j) << std::endl;
					}});
				std::thread th4([&v4, &V4, &S]() {
					if (V4.activ) {
						v4 = S.mirrorReflectivity(V4);
									//std::cout << std::to_string(res.E[i][j].T.x) << '\t' << std::to_string(res.E[i][j].T.y) << '\t' << std::to_string(res.E[i][j].T.z) << std::endl;
					}
					else {
						v4 = V4;
									//std::cout << std::to_string(i * 10 + j) << std::endl;
					}});
				th1.join();
				th2.join();
				th3.join();
				th4.join();
				res.E[i][j * 4] = v1;
				res.E[i][j * 4 + 1] = v2;
				res.E[i][j * 4 + 2] = v3;
				res.E[i][j * 4 + 3] = v4;

			}
		}

		return res;
	}
	emiter reflectToMirror(axisSpherSurface S) {
		emiter res;
		ray v1;
		ray v2;
		ray v3;
		ray v4;
		ray V1;
		ray V2;
		ray V3;
		ray V4;
		
		int reduceNumberRay = numberRay / 4;
		for (int i = 0; i < numberRay;i++) {
			for (int j = 0; j < reduceNumberRay;j++) {
				V1 = this->E[i][j * 4];
				V2 = this->E[i][j * 4 + 1];
				V3 = this->E[i][j * 4 + 2];
				V4 = this->E[i][j * 4 + 3];
				//std::cout << std::to_string(i) << '\t' << std::to_string(j) << std::endl;
				std::thread th1([&v1, &V1, &S]() {
					if (V1.activ) {
						v1 = S.mirrorReflectivity(V1);
						//std::cout << std::to_string(res.E[i][j].T.x) << '\t' << std::to_string(res.E[i][j].T.y) << '\t' << std::to_string(res.E[i][j].T.z) << std::endl;
					}
					else {
						v1 = V1;
						//std::cout << std::to_string(i * 10 + j) << std::endl;
					}});
					std::thread th2([&v2, &V2, &S]() {
						if (V2.activ) {
							v2 = S.mirrorReflectivity(V2);
							//std::cout << std::to_string(res.E[i][j].T.x) << '\t' << std::to_string(res.E[i][j].T.y) << '\t' << std::to_string(res.E[i][j].T.z) << std::endl;
						}
						else {
							v2 = V2;
							//std::cout << std::to_string(i * 10 + j) << std::endl;
						}});
						std::thread th3([&v3, &V3, &S]() {
							if (V3.activ) {
								v3 = S.mirrorReflectivity(V3);
								//std::cout << std::to_string(res.E[i][j].T.x) << '\t' << std::to_string(res.E[i][j].T.y) << '\t' << std::to_string(res.E[i][j].T.z) << std::endl;
							}
							else {
								v3 = V3;
								//std::cout << std::to_string(i * 10 + j) << std::endl;
							}});
							std::thread th4([&v4, &V4, &S]() {
								if (V4.activ) {
									v4 = S.mirrorReflectivity(V4);
									//std::cout << std::to_string(res.E[i][j].T.x) << '\t' << std::to_string(res.E[i][j].T.y) << '\t' << std::to_string(res.E[i][j].T.z) << std::endl;
								}
								else {
									v4 = V4;
									//std::cout << std::to_string(i * 10 + j) << std::endl;
								}});
								th1.join();
								th2.join();
								th3.join();
								th4.join();
								res.E[i][j * 4] = v1;
								res.E[i][j * 4 + 1] = v2;
								res.E[i][j * 4 + 2] = v3;
								res.E[i][j * 4 + 3] = v4;

			}
			
		}
		return res;
		}
	emiter translateToPlane(plane P) {
		emiter res;
		bool error;
		for (int i = 0; i < numberRay;i++) {
			for (int j = 0; j < numberRay;j++) {
				if (this->E[i][j].activ) {
					res.E[i][j].R = P.transvers(this->E[i][j], error);
					//std::cout << std::to_string(this->E[i][j].T.x) << '\t' << std::to_string(this->E[i][j].T.y) << '\t' << std::to_string(this->E[i][j].T.z) << std::endl;
					res.E[i][j].T = this->E[i][j].T;
					if (error) {
						res.E[i][j] = this->E[i][j];
						res.E[i][j].activ = false; 
					}
				}
				else {
					res.E[i][j] = this->E[i][j];
				}
			}
		}
		return res;
	}
	emiter refractiv(axisSpherSurface S, double n1, double n2) {
		emiter res;
		ray v1;
		ray v2;
		ray v3;
		ray v4;
		ray V1;
		ray V2;
		ray V3;
		ray V4;

		int reduceNumberRay = numberRay / 4;
		for (int i = 0; i < numberRay;i++) {
			for (int j = 0; j < reduceNumberRay;j++) {
				V1 = this->E[i][j * 4];
				V2 = this->E[i][j * 4 + 1];
				V3 = this->E[i][j * 4 + 2];
				V4 = this->E[i][j * 4 + 3];
				//std::cout << std::to_string(i) << '\t' << std::to_string(j) << std::endl;
				std::thread th1([&v1, &V1, &S, &n1, &n2]() {
					if (V1.activ) {
						bool PVO = false;
						v1 = S.refractive(V1,n1,n2, PVO);
						if (PVO) {
							v1.activ = false;
						}
						//std::cout << std::to_string(res.E[i][j].T.x) << '\t' << std::to_string(res.E[i][j].T.y) << '\t' << std::to_string(res.E[i][j].T.z) << std::endl;
					}
					else {
						v1 = V1;
						//std::cout << std::to_string(i * 10 + j) << std::endl;
					}});
					std::thread th2([&v2, &V2, &S, &n1, &n2]() {
						if (V2.activ) {
							bool PVO = false;
							v2 = S.refractive(V2, n1, n2, PVO);
							if (PVO) {
								v2.activ = false;
							}
							//std::cout << std::to_string(res.E[i][j].T.x) << '\t' << std::to_string(res.E[i][j].T.y) << '\t' << std::to_string(res.E[i][j].T.z) << std::endl;
						}
						else {
							v2 = V2;
							//std::cout << std::to_string(i * 10 + j) << std::endl;
						}});
						std::thread th3([&v3, &V3, &S, &n1, &n2]() {
							if (V3.activ) {
								bool PVO = false;
								v3 = S.refractive(V3, n1, n2, PVO);
								if (PVO) {
									v3.activ = false;
								}
								//std::cout << std::to_string(res.E[i][j].T.x) << '\t' << std::to_string(res.E[i][j].T.y) << '\t' << std::to_string(res.E[i][j].T.z) << std::endl;
							}
							else {
								v3 = V3;
								//std::cout << std::to_string(i * 10 + j) << std::endl;
							}});
							std::thread th4([&v4, &V4, &S, &n1, &n2]() {
								if (V4.activ) {
									bool PVO = false;
									v4 = S.refractive(V4, n1, n2, PVO);
									if (PVO) {
										v4.activ = false;
									}
									//std::cout << std::to_string(res.E[i][j].T.x) << '\t' << std::to_string(res.E[i][j].T.y) << '\t' << std::to_string(res.E[i][j].T.z) << std::endl;
								}
								else {
									v4 = V4;
									//std::cout << std::to_string(i * 10 + j) << std::endl;
								}});
								th1.join();
								th2.join();
								th3.join();
								th4.join();
								res.E[i][j * 4] = v1;
								res.E[i][j * 4 + 1] = v2;
								res.E[i][j * 4 + 2] = v3;
								res.E[i][j * 4 + 3] = v4;
			}
		}
		return res;
	}
};