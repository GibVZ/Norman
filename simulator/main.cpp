#include<iostream> // вывод в консоль
#include<cstdint>  // чтобы он не ругался на uint32_t
#include<cmath>    // математика pow и sqrt
#include<fstream>  // Вывод файл
// using namespace std;

// precision type is double or float
using prectype = double; // переименовываем double, чтобы если захотелось поменять точность не надо было заменять все double'ы
using u32 = uint32_t;    // используем unsigned потому что это для количемтва частиц


const prectype sigma = 0.1; //         m
const prectype epsilon = 1e-1; //      J

const prectype left_border  = 0; //    m
const prectype right_border = 1; //    m
const prectype Vmax = 10; //        m / s


class Vector3{          // Класс, отвечающий за вектора - чтобы каждый раз не прописывать покомпонетное сложение/умножения на скаляр
	public:
		prectype coords[3];
		
		Vector3(prectype x, prectype y, prectype z){
			coords[0] = x;
			coords[1] = y;
			coords[2] = z;
		}
		
		Vector3(){
			coords[0] = 0;
			coords[1] = 0;
			coords[2] = 0;
		}
		
		Vector3(Vector3& other){
			coords[0] = other.coords[0];
			coords[1] = other.coords[1];
			coords[2] = other.coords[2];
		}
/*		
		Vector3& operator=(Vector3& other) {
			coords[0] = other.coords[0];
			coords[1] = other.coords[1];
			coords[2] = other.coords[2];	
			return *this;
		}
*/		
		Vector3& operator=(Vector3 other) {
			coords[0] = other.coords[0];
			coords[1] = other.coords[1];
			coords[2] = other.coords[2];	
			return *this;
		}
		
		prectype& operator[](int ind){
			return coords[ind];
		}
		
		prectype norm(){
			prectype sum = 0;
			for(int i = 0; i < 3; i++){
				sum += coords[i]*coords[i];
			}
			return std::sqrt(sum);
		}
		
		prectype norm2() {
			prectype sum = 0;
			for(int i = 0; i < 3; i++){
				sum += coords[i]*coords[i];
			}
			return sum;
		}
};

Vector3 operator+(Vector3 self, Vector3 other){ // определение сложения для вышеописанного класса
	Vector3 ret(0., 0., 0.);
	for(int i = 0; i < 3; i++){
		ret[i] = self[i] + other[i];
	}
	return ret;
}

Vector3 operator-(Vector3 self, Vector3 other){ // вычитание
	Vector3 ret(0., 0., 0.);
	for(int i = 0; i < 3; i++){
		ret[i] = self[i] - other[i];
	}
	return ret;
}

Vector3 operator*(prectype scalar, Vector3 self){ // умнодение на скляр справа
	Vector3 ret(0., 0., 0.);
	for(int i = 0; i < 3; i++){
		ret[i] = self[i] * scalar;
	}
	return ret;
}

Vector3 operator*(Vector3 self, prectype scalar){ // умнодение на скаляр слева
	Vector3 ret(0., 0., 0.);
	for(int i = 0; i < 3; i++){
		ret[i] = self[i] * scalar;
	}
	return ret;
}

Vector3 rand_pos(){     // функция для начальных положений частиц(на случай, если захочется поменять распределение или еще что-нибудь)
	Vector3 ret;
	for(int i = 0; i < 3; i++){
		ret[i] = 2 * (prectype)std::rand()/RAND_MAX - 1;
	}
	return ret;
}

Vector3 rand_vel(){     // аналогично со скоростью
	Vector3 ret;
	for(int i = 0; i < 3; i++){
		ret[i] = 2 * (prectype)std::rand() / RAND_MAX - 1;
	}
	return ret  * (Vmax / ret.norm());
}

class Particle{              // Класс частицы
	public: // Уровень доступа
		Vector3 position, velocity, force; // координаты, скорость и сила соответственно
		prectype k_energy, g_energy;
		prectype mass; // масса частицы
		Particle(){ // Здесь задаются начальные условия(конструктор класса, который вызывается при создании экземпляра)
			mass = 1;
			k_energy = 0;
			g_energy = 0;
			position = rand_pos(); // случайное положение (третья координата ноль, чтобы можно было визуализировать в pyagme(2D))
			velocity = rand_vel(); // аналагично для скорости
		}
		// friend Vector3 operator+(Vector3, Vector3);
};



Vector3 Force(Particle self, Particle other){
	prectype r = (other.position - self.position).norm2();
	return 24 * epsilon * (std::pow(sigma, 6) / std::pow(r, 4) - 2 * std::pow(sigma, 12) / std::pow(r, 7)) * (other.position - self.position);  
}
/*
Vector3 Force(Particle self, Particle other){ // Функция силы взаиможействия двух частиц, если захочется ее поменять(потенциал леннарда-джонса в идеале)
	return (self.position - other.position) * (1/std::pow((self.position - other.position).norm(), 3));
}
*/

void compute_energy(size_t particles_amount, Particle *particles) {
	for (size_t i = 0; i < particles_amount; ++i) {
		particles[i].k_energy = particles[i].mass * particles[i].velocity.norm2() / 2;
		particles[i].g_energy = 0;
		for (size_t j = 0; j < i; ++j) {
			Particle& p1 = particles[i];
			Particle& p2 = particles[j];
			prectype r = (p1.position - p2.position).norm();
			prectype ge = 4 * epsilon * (std::pow(sigma / r, 12) - std::pow(sigma / r, 6));
			p1.g_energy += ge;
			p2.g_energy += ge;
		}
	}
}

void compute_forces(size_t particles_num, Particle* particles){ // вычисляем в моменте силу, действующую на каждую из частиц
	for(size_t i = 0; i < particles_num; i++){
		particles[i].force = Vector3(0, 0, 0); // Если не работает - попробовать обнулять силу в отдельном цикле
		for(size_t j = 0; j < i; j++){
			Vector3 force = Force(particles[i], particles[j]);
			particles[i].force = particles[i].force + force;
			particles[j].force = particles[j].force - force;
		}
	}
}

// обновляются скорости с учетом действующий сил
void compute_velocities(size_t particles_num, Particle* particles, prectype delta_t){
	for(size_t i = 0; i < particles_num; i++){
		particles[i].velocity = particles[i].velocity + particles[i].force * (delta_t / particles[i].mass);
	}
}


prectype position_reflect(prectype position){
	/*position = -position;
	while(std::abs(position) > 1)
		position--;
	while(position < 0)
		position++;
	return position;*/
	if(position < left_border)
		return 2 * left_border - position;
	else if(position > right_border)
		return 2 * right_border - position;
	else
		return position;
}

// обновляются координаты частиц
void compute_positions(size_t particles_num, Particle* particles, prectype delta_t){
	for(size_t i = 0; i < particles_num; i++){
		Vector3 new_position = particles[i].position + particles[i].velocity * delta_t;
		for(size_t k = 0; k < 3; k++){
			if(new_position[k] != position_reflect(new_position[k])){
				new_position[k] = position_reflect(new_position[k]); 
				particles[i].velocity[k] = -particles[i].velocity[k];
			}
		}
		particles[i].position = new_position;
	}
}

Particle* grid_init(size_t particles_amount) {
	size_t arow = std::ceil(std::pow(particles_amount, 1./3));
	// std::cout << arow << "\n";
	Particle *ret = new Particle[particles_amount];
	for (size_t i = 0; i < arow; ++i) {
		for (size_t j = 0; j < arow; ++j) {
			for (size_t k = 0; k < arow; ++k) {
				// std::cout << x << " " << y << " " << z << "\n";
				if (particles_amount <= (i * arow + j) * arow + k)
					break;
				prectype x = (right_border - left_border) * double(i + 1) / (arow + 1) + left_border;
				prectype y = (right_border - left_border) * double(j + 1) / (arow + 1) + left_border;
				prectype z = (right_border - left_border) * double(k + 1) / (arow + 1) + left_border;
				ret[(i * arow + j) * arow + k].position = Vector3(x, y, z);
			}
		}
	}
	return ret;
}

int main(int argc, char* argv[]){
	if(argc <= 4){
		std::cout << "You must give 4 arguments(output file name, delta_t, particles amount, steps amount) but you give " << argc;
		return 1;
	}
	std::srand(1);
	size_t particles_amount = std::atoi(argv[3]);
	size_t steps = std::atoi(argv[4]);
	prectype delta_t = std::atof(argv[2]);
	Particle* particles = grid_init(particles_amount);// new Particle[particles_amount];
	
	std::ofstream os(argv[1]);
	if(os.is_open()){
		for(size_t i = 0; i < particles_amount; i++){
			for(size_t j = 0; j < 3; j++)
				os << particles[i].position[j] << " ";
			//os << ";";
			for(size_t j = 0; j < 3; j++)
				os << particles[i].velocity[j] << " ";
			os << particles[i].k_energy << " " << particles[i].g_energy << "\n";
		}
		os << "\n";
	}
	
	for(size_t step = 0; step < steps; step++){
		compute_energy(particles_amount, particles);
		compute_forces(particles_amount, particles);
		compute_velocities(particles_amount, particles, delta_t);
		compute_positions(particles_amount, particles, delta_t);
		// std::cout << step << "iteration was finished\n";
		if(os.is_open()){
			for(size_t i = 0; i < particles_amount; i++){
				for(size_t j = 0; j < 3; j++)
					os << particles[i].position[j] << " ";
				//os << ";";
				for(size_t j = 0; j < 3; j++)
					os << particles[i].velocity[j] << " ";
				os << particles[i].k_energy << " " << particles[i].g_energy << "\n";
			}
			os << "\n";
		}
	}
	
	delete[] particles;
}
