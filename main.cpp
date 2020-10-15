

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <cstring>
#include <iostream>
#include <algorithm>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define EPSILON 1E-250

template<typename T> T sqr(const T x) { return x*x; };


class Gaussian2DKernel {
public:
	Gaussian2DKernel(double regularization = 25., int W = 0, int H = 0) {
		gamma = regularization;
		W_ = W;
		H_ = H;
		N = W_*H_;
		kernel1d = new double[std::max(W, H)];
		for (int i=0; i<std::max(W, H); i++) {
			kernel1d[i] = std::max(EPSILON, exp(-i*i / gamma));
		}
	}

	Gaussian2DKernel(const Gaussian2DKernel& b) {
		gamma = b.gamma;
		W_ = b.W_;
		H_ = b.H_;
		N = W_*H_;
		kernel1d = new double[std::max(W_, H_)];
		std::memcpy(kernel1d, b.kernel1d, std::max(W_, H_)*sizeof(double));
	}

	Gaussian2DKernel& operator=(const Gaussian2DKernel& b) {
		gamma = b.gamma;
		W_ = b.W_;
		H_ = b.H_;
		N = W_*H_;
		kernel1d = new double[std::max(W_, H_)];
		std::memcpy(kernel1d, b.kernel1d, std::max(W_, H_)*sizeof(double));
		return *this;
	}

	~Gaussian2DKernel() {
		delete[] kernel1d;
	}

	double operator()(int i, int j) {
		int x0 = i%W_;
		int y0 = i/W_;
		int x1 = j%W_;
		int y1 = j/W_;
		return std::max(EPSILON, exp(-(sqr(x0-x1)+sqr(y0-y1))/gamma));
	}

	void convolve(const double* u, double* result) const {
		double* tmp = new double[W_*H_]; // allocating here ; otherwise not thread-safe

		for (int i=0; i<H_; i++) {
			for (int j=0; j<W_; j++) {
				double conv = 0;
				for (int k=0; k<W_; k++) {
					conv+=kernel1d[abs(j-k)]*u[i*W_ + k];
				}
				tmp[i+j*H_] = conv;
			}
		}

		for (int j=0; j<W_; j++) {
			for (int i=0; i<H_; i++) {
				double conv = 0;
				for (int k=0; k<H_; k++) {
					conv+=kernel1d[abs(i-k)]*tmp[k + j*H_];
				}
				result[i*W_+j] = conv;
			}
		}
		delete[] tmp;
	}

	void convolveAdjoint(const double* u, double* result) const {
		convolve(u, result);
	}

	double gamma;
	int W_, H_;
	int N;
	double* kernel1d;
};


int main() {

	int W, H, C;
	
	//stbi_set_flip_vertically_on_load(true);
	unsigned char *image1 = stbi_load("evol1.bmp",
                                 &W,
                                 &H,
                                 &C,
                                 1);
	std::vector<double> image_double1(W*H);
	for (int i=0; i<W*H; i++)
		image_double1[i] = image1[i];

	int sum1 = 0;
	for (int i=0; i<W*H; i++)
		sum1 += image_double1[i];
	for (int i=0; i<W*H; i++)
		image_double1[i] = image_double1[i]/sum1;

	std::cout << "sum1 " << sum1 << std::endl;

	unsigned char *image2 = stbi_load("evol2.bmp",
                                 &W,
                                 &H,
                                 &C,
                                 1);
	std::vector<double> image_double2(W*H);
	for (int i=0; i<W*H; i++)
		image_double2[i] = image2[i];
	
	int sum2 = 0;
	for (int i=0; i<W*H; i++)
		sum2 += image_double2[i];
	for (int i=0; i<W*H; i++)
		image_double2[i] = image_double2[i]/sum2;

	std::cout << "sum2 " << sum2 << std::endl;

	unsigned char *image3 = stbi_load("evol3.bmp",
                                 &W,
                                 &H,
                                 &C,
                                 1);
	std::vector<double> image_double3(W*H);
	for (int i=0; i<W*H; i++)
		image_double3[i] = image3[i];

	int sum3 = 0;
	for (int i=0; i<W*H; i++)
		sum3 += image_double3[i];
	for (int i=0; i<W*H; i++)
		image_double3[i] = image_double3[i]/sum3;

	std::cout << "sum3 " << sum3 << std::endl;

	std::vector<double> image_b1(W*H);
	for (int i=0; i<W*H; i++)
		image_b1[i] = 1;

	std::vector<double> image_b2(W*H);
	for (int i=0; i<W*H; i++)
		image_b2[i] = 1;

	std::vector<double> image_b3(W*H);
	for (int i=0; i<W*H; i++)
		image_b3[i] = 1;

	std::vector<double> image_a1(W*H);
	std::vector<double> image_a2(W*H);
	std::vector<double> image_a3(W*H);

	std::vector<double> P_lambda(W*H);

	std::vector<double> lamda(3);
	lamda[0]=1.0/3;
	lamda[1]=1.0/3;
	lamda[2]=1.0/3;

	Gaussian2DKernel K(25,W,H);

	int L = 20;
	for (int k = 0; k< L; k++){
		

		K.convolve(&image_b1[0],&image_b1[0]);
		K.convolve(&image_b2[0],&image_b2[0]);
		K.convolve(&image_b3[0],&image_b3[0]);

		for (int i=0; i<W*H; i++){
			image_a1[i]= image_double1[i]/(std::max(EPSILON, image_b1[i]));
			image_a2[i]= image_double2[i]/(std::max(EPSILON, image_b2[i]));
			image_a3[i]= image_double3[i]/(std::max(EPSILON, image_b3[i]));
		}

		K.convolve(&image_a1[0],&image_a1[0]);
		K.convolve(&image_a2[0],&image_a2[0]);
		K.convolve(&image_a3[0],&image_a3[0]);

		for (int i=0; i<W*H; i++){
			P_lambda[i] = pow(image_a1[i],lamda[0]);
			P_lambda[i] = P_lambda[i]*pow(image_a2[i],lamda[1]);
			P_lambda[i] = P_lambda[i]*pow(image_a3[i],lamda[2]);
		}

		// K.convolve(&image_a1[0],&image_a1[0]);
		// K.convolve(&image_a2[0],&image_a2[0]);
		// K.convolve(&image_a3[0],&image_a3[0]);

		for (int i=0; i<W*H; i++){
			image_b1[i]= P_lambda[i]/(std::max(EPSILON, image_a1[i]));
			image_b2[i]= P_lambda[i]/(std::max(EPSILON, image_a2[i]));
			image_b3[i]= P_lambda[i]/(std::max(EPSILON, image_a3[i]));
		}

	}

	for (int i=0; i<W*H; i++){
		P_lambda[i] = P_lambda[i] * (sum1 * lamda[0] + sum2 * lamda[1] + sum3 * lamda[2]);
		if (P_lambda[i] != 0){
			std::cout << "lamda " << P_lambda[i]<< std::endl;
		}
	}
	std::vector<unsigned char> P_lambda_result(W*H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			P_lambda_result[i*W + j] = P_lambda[i*W+j];
		}
	}
	stbi_write_png("Lambda_result.png", W, H, 1, &P_lambda_result[0], 0);




	return 0;
}