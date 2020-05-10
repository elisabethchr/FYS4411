#include "stochasticgradientdescent.h"
#include <cassert>
#include <iostream>
#include <vector>
#include <armadillo>
#include "optimizer.h"
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"
#include "../Hamiltonians/hamiltonian.h"


StochasticGradientDescent::StochasticGradientDescent(System* system, double eta) :
    Optimizer(system){
    m_eta = eta;
}

/* Calculate optimized weights via stochastic gradient descent */
void StochasticGradientDescent::computeWeights(arma::vec gradE){
    int m_nv = m_system->getNumberVisibleNodes();
    int m_nh = m_system->getNumberHiddenNodes();

    m_a = m_system->getWaveFunction()->get_a();
    m_b = m_system->getWaveFunction()->get_b();
    m_w = m_system->getWaveFunction()->get_w();

    for (int i=0; i<m_nv; i++){
        m_a[i] = m_a[i] - m_eta*gradE[i];
    }

    for (int j=0; j<m_nh; j++){
        m_b[j] = m_b[j] - m_eta*gradE[m_nv + j];
    }

    int count = m_nh + m_nv;
    for (int i=0; i<m_nv; i++){
        for (int j=0; j<m_nh; j++){
            m_w(i, j) = m_w(i, j) - m_eta*gradE[count];
            count++;
        }
    }

    m_system->getWaveFunction()->set_a(m_a);
    m_system->getWaveFunction()->set_b(m_b);
    m_system->getWaveFunction()->set_w(m_w);
}
