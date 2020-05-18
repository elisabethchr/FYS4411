#pragma once
#include <vector>
#include <armadillo>

class WaveFunction {
public:
    WaveFunction(class System* system);

    virtual arma::vec set_X(arma::vec X) = 0;    // set visible nodes vector (i.e. position vector)
    virtual arma::vec set_a(arma::vec a) = 0;    // set visible bias vector
    virtual arma::vec set_b(arma::vec b) = 0;    // set hidden bias vector
    virtual arma::mat set_w(arma::mat w) = 0;    // set interaction of biases matrix

    virtual arma::vec get_X() = 0;    // get visible nodes vector (i.e. position vector)
    virtual arma::vec get_a() = 0;    // get visible bias vector
    virtual arma::vec get_b() = 0;    // get hidden bias vector
    virtual arma::mat get_w() = 0;    // get interaction of biases matrix

    virtual double evaluate(arma::vec position) = 0;

    virtual double computeDerivative_analytic(arma::vec position, int coor) = 0;
    virtual double computeDoubleDerivative_analytic() = 0;
    virtual double getDistance(int p, int q) = 0;
    virtual double getSigma() = 0;
    virtual double getInitializationInterval() = 0;
    virtual bool getGaussian() = 0;

protected:
    class System* m_system = nullptr;

    int m_nv = 0;       // number visible nodes
    int m_nh = 0;       // number hidden nodes
    int m_part = 0;     // number of particles
    int m_dim = 0;      // number of dimensions
    std::mt19937_64 m_randomEngine; // for the distributions

private:

    virtual void setupInitialState() = 0;

    virtual double sigmoid(double x) = 0; // the logistic function
    virtual double v(int j) = 0;          // for simplification
    virtual arma::mat hadamardProd(arma::mat w) = 0; // calculates the hadamard product of w with itself

};
