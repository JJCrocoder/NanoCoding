//We define the random generator functions
float uniform(float min, float max) { // Random real number in an uniform distribution
    return min + (max-min)*rand()/RAND_MAX;
}

int rand_num(int min, int max) { // Random integer number in an uniform discrete distribution
    return rand() % (max-min+1) + min;
}

float randn(float mean, float var) {// Random real number in a gausian distribution
	default_random_engine generator; // we define the number generator
	normal_distribution<double> distribution(mean,var); // We define the distribution
	return distribution(generator);
}
