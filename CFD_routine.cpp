double waveform_baseline(TH1D* waveform) {
    double sample_sum = 0.;
    int count = 0;
    for (int k=1; k <20; k++) {
        sample_sum += waveform->GetBinContent(k);
        count++;
    }
    double baseline = sample_sum/float(count);
    return baseline;
}

int global_timing(TH1D* waveform, const int number_samples, double baseline) {
        // Find global minimum
        double ymin = 99999999.;
        int imin = 0;
        for (int k = 1; k < number_samples-2; k++){
            if (waveform->GetBinContent(k)-baseline < ymin) {
                ymin = waveform->GetBinContent(k)-baseline;
                imin = k;
                // std::cout<< "   " << ymin << std::endl;
            }
        }
    return imin;
}


float CFD_timing(TH1D* waveform, double baseline, const int number_samples, const int global_imin, const float startp, const float endp, const float percentage) {
    double y_min = waveform->GetBinContent(global_imin)-baseline;
    double y_end = endp*y_min;
    double y_start = startp*y_min;
    double rise_amplitude = (endp-startp)*y_min;
    int j = global_imin;
    int j_end = 0;
    while (((waveform->GetBinContent(j)-baseline) < y_end) && (j > 1)) {
        j--;
        j_end = j-1;
    }
    j = global_imin;
    int j_start = 0;
    while (((waveform->GetBinContent(j)-baseline) < y_start) && (j > 1)) {
        j--;
        j_start = j+1;
    }
    double b = (y_end-y_start)/(double(j_end)-double(j_start));
    double a = y_start - j_start*b;
    double iCFD = (percentage*rise_amplitude-a)/b;
    return iCFD;
}

double charge(TH1D* waveform, double baseline, const int A, const int B, double verticalScaleFactor, double horizontalScaleFactor) {
    double charge = 0.;
    for (int j = A; j < B; j++) {
        charge += (waveform->GetBinContent(j)-baseline)*verticalScaleFactor*horizontalScaleFactor;
    }
    return charge;
}

// Conversion factor
#define mppc_DIGITIZER_FULL_SCALE_RANGE 2.5 // Vpp
#define mppc_DIGITIZER_RESOLUTION       12 // bits
#define mppc_DIGITIZER_SAMPLE_RATE      80000000 // 
double digiCounts = pow(2.0, mppc_DIGITIZER_RESOLUTION); // counts
double verticalScaleFactor = 1.e3*mppc_DIGITIZER_FULL_SCALE_RANGE/digiCounts; // V / counts
double horizontalScaleFactor = 1.e9/mppc_DIGITIZER_SAMPLE_RATE; // ns/Sa

// Cast data into TH1D* histograms
TH1D* hwaveform = (TH1D*) h_ampl[j];
double baseline = waveform_baseline(hwaveform);
int global_sample = global_timing(hwaveform, number_samples, baseline);
double famplitude = (hwaveform->GetBinContent(global_sample)-baseline)*verticalScaleFactor;
float CFD_sample  = CFD_timing(hwaveform, baseline, number_samples, global_sample, start_rise_percentage, end_rise_percentage, percentage);
double pulse_charge = charge(hwaveform, baseline, 1, number_samples, verticalScaleFactor, horizontalScaleFactor);