#include "initalizestate.h"
#include <iostream>
#include <cassert>
#include <armadillo>
#include "Math/random.h"
#include "../particle.h"
#include "../system.h"

using std::cout;
using std::endl;
using namespace std;
using namespace arma;

InitializeState::InitializeState(System* system, int n_hidden, int n_visible, int dim) : InitialState(system) {
    assert(n_hidden > 0 && n_visible > 0);
    m_nh = n_hidden;
    m_nv  = n_visible;
    m_dim = dim;

    std::random_device rd;
    m_randomEngine = std::mt19937_64(rd());

    /* The InitialState class is in charge of everything to do with the
     * initialization of the system; this includes determining the number of
     * hidden and visible nodes used. To make sure everything
     * works as intended, this information is passed to the system here.
     */

    m_system->setNumberHiddenNodes(n_hidden);
    m_system->setNumberVisibleNodes(n_visible);

    setupInitialState();

}

void InitializeState::setupInitialState(){
    std::vector<double> m_x;
    std::vector<double> m_a;
    std::vector<double> m_b;
    arma::mat m_w = mat(m_nh, m_nv);

    std::uniform_real_distribution<double> uniform(-1.0,1.0);

    for (int i=0; i<m_nv; i++){
        m_x.push_back(uniform(m_randomEngine));
        m_a.push_back(uniform(m_randomEngine));
        m_b.push_back(uniform(m_randomEngine));
        for (int j=0; j<m_nh; j++){
            m_w(i, j) = uniform(m_randomEngine);
        }
    }

    setupPositions();
}

//void InitializeState::setupWeights(){}

void InitializeState::setupPositions(){
    std::uniform_real_distribution<double> uniform_initX(-0.5,0.5);

    std::vector<double> m_x;
    for (int i=0; i<m_nv; i++){
        m_x.push_back(uniform_initX(m_randomEngine));
    }
}
