#include "CorrectStepChanges.h"

CorrectStepChanges::CorrectStepChanges():Tool(){}


bool CorrectStepChanges::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_variables.Get("verbosity", m_verbose);
	m_verbose = 10;
	
	// All F-distribution probabilities above this value will pass the
	// variance consistency test in ze3ra_baseline(). That is, p_critical
	// is the maximum p-value for which we will reject the null hypothesis
	// of equal errs.
	m_variables.Get("PCritical", p_critical);
	
	m_variables.Get("stddev_tolerance",stddev_tolerance);
	m_variables.Get("conc_change_tolerance",conc_change_tolerance);
	m_variables.Get("sample_size",sample_size);
	m_variables.Get("ignored_pts",ignored_pts);
	
	m_data= &data;
	
	return true;
}

// FIXME perhaps an alternative to adding a constant offset to bring the new concentration values in line
// would be to update the absorption region shape - particularly if the fit chi2 is worse now than before?

bool CorrectStepChanges::Execute(){
	
	// see if we have new data to add to DB
	std::string ledname="";
	get_ok = m_data->CStore.Get("NewMatthewAnalyse",ledname);
	
	// clear previous results if applicable
	m_data->CStore.Remove("NewCorrectedConc");
	
	// do we have a new measurement?
	if(!get_ok || ledname==""){
		Log("SaveToDB::CorrectStepChanges no new measurement",v_debug,verbosity);
		
		// this shouldn't get called in the AnalyseOldFiles ToolChain because
		// there's a new measurement on every loop. But for data,
		// this will be called on the loop after a measurement has been taken.
		// on the first such loop after a measurement, retroactively grab
		//  the measurement number and timestamp (generated/used by SaveToDB)
		if(new_meas){
			get_ok = m_data->CStore.Get("last_measurement_num",thismeasurementnum);
			if(!get_ok){
				// can't think why this would happen except if the SaveToDB Tool
				// isn't in the toolchain (in which case it doesn't matter) 
				// or crashed out (in which case the original measurement wasn't saved?)
				// in any case, we should try to guess, or maintain monotonicity...?
				++thismeasurementnum;
			}
			measurementnum.push_back(thismeasurementnum);
			
			// same for the timestamp
			get_ok = m_data->CStore.Get("last_measurement_timestamp",thistimestamp);
			if(!get_ok){
				// try to fall back to something that will retain monotonicity of timestamps
				thistimestamp=to_simple_string(m_data->measurment_time);
			}
			timestamps.push_back(thistimestamp);
			
			// ok, got the values we needed.
			new_meas = false;
		}
		
		return true;
	}
	// else a new concentration measurement has been taken. We need to follow suit.
	new_meas = true;
	
	// shuffle values along
	// 'corrected_set's store both conc and original error
	if(corrected_set2.size()==sample_size){
		corrected_set1.push_back(corrected_set2.front());
		corrected_set2.pop_front();
	}
	
	if(corrected_set1.size()>(sample_size+ignored_pts)){
		// a datapoint has made it through the pipeline and is ready for writing out.
		// put corrected value into CStore
		m_data->CStore.Set("NewCorrectedConc", lednames.front());
		lednames.pop_front();
		m_data->CStore.Set("conc_and_err_corr", corrected_set1.front());
		corrected_set1.pop_front();
		m_data->CStore.Set("corr_conc_meas_num",measurementnum.front()); // for matching with original
		measurementnum.pop_front();
		m_data->CStore.Set("corr_conc_timestamp", timestamps.front());
		timestamps.pop_front();
		m_data->CStore.Set("step_vetoed",step_vetoed.front());
		step_vetoed.pop_front();
		m_data->CStore.Set("step_found",step_detected.front());
		step_detected.pop_front();
		m_data->CStore.Set("step_applied",step_applied.front());
		step_applied.pop_front();
		m_data->CStore.Set("accum_step_changes", running_step.front());
		running_step.pop_front();
	}
	
	// 'gdconcs_set's store only concentration for calculation of mean & stddev
	if(gdconcs_set2.size()==sample_size){
		gdconcs_set1.push_back(gdconcs_set2.front());
		gdconcs_set2.pop_front();
	}
	
	if(gdconcs_set1.size()>(sample_size+ignored_pts)){
		gdconcs_set1.pop_front();
	}
	
	// get the new datapoint
	std::pair<double,double> conc_and_err;
	get_ok = m_data->CStore.Get("conc_and_err",conc_and_err);
	
	// if we didn't get one...
	if(!get_ok){
		Log(m_unique_name+" Error! No conc_and_err in CStore!",v_error,m_verbose);
		
		// we always need to maintain one datapoint per execute to stay in sync
		corrected_set2.push_back(std::pair<double,double>{0,0});
		
		// ~consider that this invalidates our step change procedure; start over~
		// probably unnecessary, the fit chi2 will probably be bad enough with a spurious datapoint anyway.
		//gdconcs_set1.clear();
		//gdconcs_set2.clear();
		gdconcs_set1.push_back(std::pair<double,double>{0,0});
		gdconcs_set2.push_back(conc_and_err.first);
		
		return true;  // not an error on our side
	}
	
	// otherwise start off by adding our corrections accumulated so far to it
	// to bring it in line with our buffered values
	conc_and_err.first += running_change;
	
	// and add it to our buffers
	corrected_set2.push_back(conc_and_err);
	gdconcs_set2.push_back(conc_and_err.first);
	
	// if we don't have enough good points to do a step check yet, just return nothing
	if(gdconcs_set1.size()<(sample_size+ignored_pts)){
		Log(m_unique_name+ " npoints: "+std::to_string(gdconcs_set1.size())
		   +": not enough",v_debug,m_verbose);
		return true;
	}
	
	// debug prints
	/*
	std::cout<<"set 1: {";
	for(int i=0; i<sample_size; ++i){
		if(i>0) std::cout<<", ";
		std::cout<<gdconcs_set1.at(i);
	}
	std::cout<<"}"<<std::endl
	         <<" ignored points: {";
	for(int i=sample_size; i<(sample_size+ignored_pts); ++i){
		if(i>sample_size) std::cout<<", ";
		std::cout<<gdconcs_set1.at(i);
	}
	std::cout<<"}\nsample 2: {";
	for(int i=0; i<sample_size; ++i){
		if(i>0) std::cout<<", ";
		std::cout<<gdconcs_set2.at(i);
	}
	std::cout<<"}"<<std::endl;
	*/
	
	// calculate mean and variance of each dataset
	double mean2 = TMath::GeomMean(gdconcs_set2.begin(),gdconcs_set2.end());
	double mean1 = TMath::GeomMean(gdconcs_set1.begin(),gdconcs_set1.end()-ignored_pts);
	double stddev2 = TMath::StdDev(gdconcs_set2.begin(),gdconcs_set2.end());
	double stddev1 = TMath::StdDev(gdconcs_set1.begin(),gdconcs_set1.end()-ignored_pts);
	
	Log(m_unique_name+" mean1: "+std::to_string(mean1)+", mean2: "+std::to_string(mean2)
	    +", stddev1: "+std::to_string(stddev1)+", stddev2: "+std::to_string(stddev2),
	    v_debug,m_verbose);
	
	// check that data is stable within both regions
	bool veto_correction = false;
	if(stddev1>stddev_tolerance){
		Log(m_unique_name+" stddev of data set 1 ("+std::to_string(stddev1)
		    +") greater than tolerance of "+std::to_string(stddev_tolerance)
		    +"; step removal not applied",v_warning,m_verbose);
		veto_correction = true;
	}
	if(stddev2>stddev_tolerance){
		Log(m_unique_name+" stddev of data set 2 ("+std::to_string(stddev2)
		    +") greater than tolerance of "+std::to_string(stddev_tolerance)
		    +"; step removal not applied",v_warning,m_verbose);
		veto_correction = true;
	}
	
	m_data->CStore.Set("StepVetoed",veto_correction);
	
	// FIXME could also track the chi2 of the fits used for these measurements
	// (or perhaps the errors on them) and only take action when the measurements are good
	
	// FIXME could also somehow take into account when change between two LEDs is not the same?
	// more likely to be false in that case...but does depend on linearity between both LEDs being equal.
	
	double this_step = 0;
	bool this_step_found = false;
	
	// check whether both datasets are consistent with each other
	Log(m_unique_name+" diff between datasets: "+std::to_string(mean2-mean1),
	    v_debug,m_verbose);
	if((std::abs(mean2-mean1) > conc_change_tolerance)){
		
		// Note that we saw a step
		this_step_found = true;
		m_data->CStore.Set("StepDetected",this_step_found);
		
		// do not apply correction if data is unstable
		if(veto_correction) return true;
		
		// otherwise remove it
		this_step = mean1 - mean2;
		m_data->CStore.Set("StepChange",this_step);
		Log(m_unique_name+" step correction this execute: "+std::to_string(this_step),v_debug,m_verbose);
		
		running_change += this_step;
		m_data->CStore.Set("SumStepChange", running_change);
		Log(m_unique_name+" accumulated steps so far: "+std::to_string(running_change),v_debug,m_verbose);
		
		// we need to retroactively fix all the concentrations of set 2
		for(int i=0; i<sample_size; ++i){
			*std::next(gdconcs_set2.begin(),i) += this_step;
			std::next(corrected_set2.begin(),i)->first += this_step;
		}
		
		// fix the ignored region. Hard to tell what to do with this,
		// as it's by definition the region between two stable regions,
		// but itself may not be stable if the 'step' is spread across multiple points.
		// as the crudest hack, just replace it with the average of the points before and after the step
		// FIXME this is not even 'doctoring' data, it's outright making it up.
		// not even going to handle the case of multiple ignored points; if we do this,
		// the fact that we have multiple identical points should highlight that we need to
		// return to handle this properly
		for(int i=sample_size; i<(sample_size+ignored_pts); ++i){
			gdconcs_set1.at(i) = 0.5*(gdconcs_set1.at(sample_size-1) + gdconcs_set2.front());
		}
		
	}
	
	lednames.push_back(ledname);
	step_vetoed.push_back(veto_correction);
	step_detected.push_back(this_step_found);
	step_applied.push_back(this_step);
	running_step.push_back(running_change);
	
	// get measurement number. if retro-actively re-analysing with the AnalyseOldFiles ToolChain,
	// we will be saving to the database on every Execute, and the LoadOldFiles Tool puts a
	// 'dbmeasurementnum' in the CStore that we can grab.
	// During normal data-taking the measurement number is generated on the fly by the SaveToDB Tool,
	// so we have to but retroactively grab it on the subsequent Execute loop (done at the top).
	get_ok = m_data->CStore.Get("dbmeasurementnum",thismeasurementnum);
	if(get_ok) measurementnum.push_back(thismeasurementnum);
	// similar situation with timestamps
	get_ok = m_data->CStore.Get("dbtimestamp",thistimestamp);
	if(get_ok) timestamps.push_back(thistimestamp);
	
	return true;
}


bool CorrectStepChanges::Finalise(){
	
	return true;
}



/*
void ADCCalibrator::ze3ra_baseline(double& baseline, double& sigma_baseline, size_t n_prev_meas){
	
	// vals, errs, and F-distribution probability values
	// ("P") for the last n_prev_meas gd concentration measurements
	std::vector<double> vals;
	std::vector<double> errs;
	
	// TODO retrieve these from database
	// SELECT values->'values'->[0] from data WHERE name='gdconcmeasure'
	// AND ledname'"+ledname+"' ORDER BY timestamp LIMIT "+std::to_string(n_prev_meas);
	vals.push_back(mean);
	errs.push_back(var);
	
	// Compute probabilities for the F-distribution test for each measurement
	std::vector<double> Ps;
	for (size_t j = 0; j < errs.size() - 1; ++j) {
		double sigma2_j = errs.at(j);
		double sigma2_jp1 = errs.at(j + 1);
		double F;
		if (sigma2_j > sigma2_jp1) F = sigma2_j / sigma2_jp1;
		else F = sigma2_jp1 / sigma2_j;
		
		double nu = (n_prev_meas - 1) / 2.;
		double P = annie_math::Regularized_Beta_Function(1. / (1. + F), nu, nu);
		
		// Two-tailed hypothesis test (we need to exclude unusually small values
		// as well as unusually large ones). The tails have equal sizes, so we
		// may use symmetry and simply multiply our earlier result by 2.
		P *= 2.;
		
		// I've never seen this problem (the numerical values for the regularized
		// beta function that I've checked all fall within [0,1]), but the book
		// Numerical Recipes includes this check in a similar block of code,
		// so I'll add it just in case.
		if (P > 1.) P = 2. - P;
		
		Ps.push_back(P);
	}
	
	// Filter out gd concentration measurements
	// whose F-distribution probability falls above the critical value
	baseline = 0.;
	double variance_baseline = 0.;
	size_t num_passing = 0;
	for (size_t k = 0; k < Ps.size(); ++k) {
		if (Ps.at(k) > p_critical) {
			++num_passing;
			baseline += vals.at(k);
			variance_baseline += errs.at(k);
		}
	}
	
	sigma_baseline = 0.;
	if (num_passing > 1) {
		baseline /= num_passing;
		
		variance_baseline *= static_cast<double>(n_prev_meas - 1)
			/ (num_passing*n_prev_meas - 1);
		// Now that we've combined the sample errs correctly, take the
		// square root to get the standard deviation
		sigma_baseline = std::sqrt( variance_baseline );
	}
	else if (num_passing == 1) {
		// We only have one sample, so all we need to
		// do is take the square root of the variance to get the standard
		// deviation.
		sigma_baseline = std::sqrt( variance_baseline );
	}
	else {
		// If no measurements passed the F-distribution test,
		// choose the one closest to passing (i.e., the one with the largest
		// P-value). For a sufficiently large number of measurements
		// (e.g., 40), such a situation should be very rare.
		// TODO: consider changing this approach
		auto max_iter = std::max_element(Ps.cbegin(), Ps.cend());
		int max_index = std::distance(Ps.cbegin(), max_iter);
		
		baseline = vals.at(max_index);
		sigma_baseline = std::sqrt( errs.at(max_index) );
	}
	
	if(m_verbose >= (v_debug+1)) {
		for ( size_t x = 0; x < Ps.size(); ++x ) {
			Log(m_unique_name+" measurement " + std::to_string(x) + ", value = "
				 + std::to_string(vals.at(x)) + ", error = "
				 + std::to_string(errs.at(x)) + ", p-value = "
				 + std::to_string(Ps.at(x)), v_debug+1, m_verbose);
		}
	}
	
	Log(m_unique_name+" found "+std::to_string(num_passing) + " measurements "
		"passing the F-test", v_debug, m_verbose);
	Log("GD conc estimate: " + std::to_string(baseline) + " ± "
		+ std::to_string(sigma_baseline) + " %", v_debug, m_verbose);
	
	return;
}
*/

/*
std::vector< CalibratedADCWaveform<double> > ADCCalibrator::make_calibrated_waveforms(const std::vector< Waveform<unsigned short> >& raw_waveforms){
	size_t n_prev_meas;
	m_variables.Get("NumBaselineSamples", n_prev_meas);
	
	// Determine the baseline for the set of raw waveforms (assumed to all
	// come from the same readout for the same channel)
	double baseline, sigma_baseline;
	ze3ra_baseline(raw_waveforms, baseline, sigma_baseline,
		n_prev_meas);
	
	std::vector< CalibratedADCWaveform<double> > calibrated_waveforms;
	for (const auto& raw_waveform : raw_waveforms) {
		
		std::vector<double> cal_data;
		const std::vector<unsigned short>& raw_data = raw_waveform.Samples();
		
		for (const auto& sample : raw_data) {
			cal_data.push_back((static_cast<double>(sample) - baseline)
				 * ADC_TO_VOLT);
		}
		
		calibrated_waveforms.emplace_back(raw_waveform.GetStartTime(),
			cal_data, baseline, sigma_baseline);
	}
	
	return calibrated_waveforms;
}
*/

