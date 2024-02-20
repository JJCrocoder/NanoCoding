


// RDF construction

		//We are not considering the last particle (it is cosidered in all steps before)
		if (i == Npart-1) continue; 
		for (int j = i+1; j < Npart; ++j) // We count the j particles for pairing
		{
			float ri[dim];	// We initialize the particle i position
			float rj[dim];	// We initialize the particle j position
			for (int k = 0; k<dim; ++k)
			{
				// Entering values from the complete list
				ri[k] = Position[i * dim + k]; 
				rj[k] = Position[j * dim + k];
			}
			dist = mic_distance(ri,rj); 	// We use the functio to calculate the distance
			if (dist>=dmax) continue;	// If the distance is larger than the limit consider we don't count
			
		}
