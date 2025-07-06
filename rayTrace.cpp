#include "ray.cpp"
const int numberOfSurface = 3;
namespace rayTrace {
	class system {
	public:
		
		double d[numberOfSurface - 1];
		emiter E[numberOfSurface];
		double n[numberOfSurface - 1];
		
		//AxisSurface surface2;
		plane entrancePupilPlane;
		plane imagePlane;
		void setEntrancePupilPlane(vector position,vector a ,vector b ) {
			this->entrancePupilPlane.set(a, b, position);
		}
		void setimagePlane(vector position, vector a, vector b) {
			this->imagePlane.set(a, b, position);
		}
		axisSpherSurface Surface[numberOfSurface - 2];
		void setSurface(double R, double d, int number) {
			//std::cout << "sdfgfg" << std::endl;
			if (number > numberOfSurface-1) {
				//std::cout << "sdfgfg" << std::endl;
				return;
			}
			/*double reald = 0;
			this->d[number - 1] = d;
			for (int i = 0;i < number;i++) {
				reald += this->d[i];
			}*/
			Surface[number].set(R, d);
		}

		/*AxisSurface surface1;
		void set1Surface(double z, double* t, double Rmax, double Rmin) {
			if (Rmax < Rmin) std::cout << "максимальный радиус должен быть больше минимального" << std::endl;
			this->surface1.setGrid(Rmin, Rmax);
			this->surface1.setZ(z);
			for (int i = 0;i < Num;i++) {
				this->surface1.t[i] = t[i];
			}
		}*/
		/*void set2Surface(double* t, double Rmax, double Rmin) {
				this->surface2.setGrid(Rmin, Rmax);
				for (int i = 0;i < Num;i++) {
					this->surface2.t[i] = t[i];
				}
		}*/
		void setEntranceEmiter(emiter e) {
			this->E[0].set(e);
		}
		ray* rayTrace(ray L) {
			if (!L.activ) { return &L; }
			static ray res[numberOfSurface - 1];
			res[0] = L;
			bool err = false;
			for (int i = 0; i < numberOfSurface - 2;i++) {
				if (err) {
					res[i + 1] = res[i];
					res[i + 1].activ = false;
				}
				bool PVO;
				res[i + 1] = this->Surface[i].refractive(res[i], this->n[i], this->n[i + 1], PVO);
				
				//std::cout << "dsfsdf" << std::endl;
			}
			return res;
			/*ray newL = surface1.mirrorReflectivity(L);
			return newL;*/
		}
		void emiterTrace() {

			for (int i = 0;i < numberOfSurface - 2;i++) {
				//std::cout << std::to_string(i);
				this->E[i + 1] = this->E[i].refractiv(this->Surface[i],this->n[i], this->n[i+1]);
				
			}
			this->E[numberOfSurface - 1] = this->E[numberOfSurface - 2].translateToPlane(imagePlane);
			//this->E[1].set(E[0].reflectToMirror(surface1));
			//this->E[2].set(this->E[1].translateToPlane(this->imagePlane));
			//this->E[2].set(E[1].reflectToMirror(surface2));
			//this->E[3].translateToPlane(this->imagePlane);

		}




	};
}