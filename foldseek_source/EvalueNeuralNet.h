//
// Created by Martin Steinegger on 19/01/2022.
//

#ifndef FOLDSEEK_EVALUENEURALNET_H
#define FOLDSEEK_EVALUENEURALNET_H
#include <cmath>
#include "kerasify/keras_model.h"
#include "BaseMatrix.h"
#include <iostream>
class EvalueNeuralNet {
private:
    BaseMatrix *subMat;
    double logDbResidueCount;
    KerasModel encoder;
    Tensor in;
    Tensor out;
public:

    EvalueNeuralNet(size_t dbResCount, BaseMatrix* subMat);

    std::pair<double, double> predictMuLambda(unsigned char * seq, unsigned int L);

    double computePvalue(double score, double lambda_, double mu) {
	double pvalue = -999;
        double h = lambda_ * (score - mu);
        if(h > 10) {
            pvalue = -h;
        } else if (h < -2.5) {
            pvalue = -exp(-exp(-h));
        } else {
            pvalue = log((1.0 - exp(-exp(-h))));
        }
	printf("@RCE@\tscore=%.3g\tlambda_=%.3g\tmu=%.3g\tpvalue=%.3g\n", score, lambda_, mu, pvalue);
	return pvalue;
    }

    double computeEvalue(double score, double lambda_, double mu){
        return exp(computePvalue(score, lambda_, mu) + logDbResidueCount);
    }
    
    double computeEvalueCorr(double score, double lambda_, double mu){
	double logPVal = computePvalue(score, lambda_, mu);    
	double dbSizeTimesLogPVal = logPVal + logDbResidueCount;
	double evalue = exp(dbSizeTimesLogPVal);
	double corrEvalue = pow(evalue, 0.32);
	printf("@RCE@\tscore=%.3g\tlambda_=%.3g\tmu=%.3g\tlogPVal=%.3g\tevalue=%.3g\tcorrEvalue=%.3g\n",
			score, lambda_, mu, logPVal, evalue, corrEvalue);
        return corrEvalue;
    }
};


#endif //FOLDSEEK_EVALUENEURALNET_H
