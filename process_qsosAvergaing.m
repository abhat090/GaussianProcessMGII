% process_qsos: run DLA detection algorithm on specified objects

% load C4 catalog
Full_catalog = ...
    load(sprintf('%s/catalog', processed_directory(release)));
% train_ind -> those LOSs without any problem (filter_flag==0) that does not have
% Civ and useful for training null model (a model without Civ absorption line)
% prior_ind -> those LOSs with Civ absorption and in the half part of test
% test_ind -> second half without any filter flag for testing null and absorption model
% on the LOSs that we know have Civ or not. So, we asses our algorithm in this way.

if (ischar(prior_ind))
    prior_ind = eval(prior_ind);
end

% My prior_ind here is already those OK sight of lines that have CIV
prior.z_qsos  = Full_catalog.all_zqso(prior_ind);
prior.c4_ind = prior_ind;
prior.z_c4 = Full_catalog.all_z_civ(prior_ind);

% filter out CIVs from prior catalog corresponding to region of spectrum below
% Ly-alpha QSO rest. In Roman's code, for detecting DLAs, instead of
% Ly-alpha they have Lyman-limit
for i = size(prior.z_c4)
    if (observed_wavelengths(civ_1548_wavelength , prior.z_c4(i)) < ...
            observed_wavelengths(min_lambda, prior.z_qsos(i)))
        prior.c4_ind(i) = false;
    end
end
prior = rmfield(prior, 'z_c4');

% enable processing specific QSOs via setting to_test_ind
if (ischar(test_ind))
    test_ind = eval(test_ind);
end

all_wavelengths    =    all_wavelengths(test_ind);
all_flux           =           all_flux(test_ind);
all_noise_variance = all_noise_variance(test_ind);
all_pixel_mask     =     all_pixel_mask(test_ind);
all_sigma_pixel    =    all_sigma_pixel(test_ind);
z_qsos             =           all_zqso(test_ind);
num_quasars        =                numel(z_qsos);

% preprocess model interpolants
% griddedInterpolant does an interpolation and gives a function handle
% based on the grided data. If the data is 2xD, {x,y} are the the same size as
% row  and columns. M is like M=f(x,y) and is like a matrix with each element
% M(i,j) = f(x(i), y(j))
mu_interpolator = ...
    griddedInterpolant(rest_wavelengths,        mu,        'linear');
M_interpolator = ...
    griddedInterpolant({rest_wavelengths, 1:k}, M,         'linear');

% initialize results with nan
min_z_c4s                   = nan(num_quasars, 1);
max_z_c4s                   = nan(num_quasars, 1);
log_priors_no_c4            = nan(num_quasars, 1);
log_priors_c4               = nan(num_quasars, 1);
log_likelihoods_no_c4       = nan(num_quasars, 1);
% sample_log_likelihoods_c4L1 = nan(num_quasars, num_C4_samples);
sample_log_likelihoods_c4L2 = nan(num_quasars, num_C4_samples);
log_likelihoods_c4L1        = nan(num_quasars, 1);
log_likelihoods_c4L2        = nan(num_quasars, 1);
log_posteriors_no_c4        = nan(num_quasars, 1);
log_posteriors_c4L1         = nan(num_quasars, 1);
log_posteriors_c4L2         = nan(num_quasars, 1);
map_N_c4L1                  = nan(num_quasars, 1);
map_N_c4L2                  = nan(num_quasars, 1);
map_z_c4L1                  = nan(num_quasars, 1);
map_z_c4L2                  = nan(num_quasars, 1);
map_sigma_c4L1              = nan(num_quasars, 1);
map_sigma_c4L2              = nan(num_quasars, 1);
sample_sigma_c4 = min_sigma + (max_sigma-min_sigma)*offset_sigma_samples;
% sample_sigma_c4 = sigma_samples;
fAVG= 1/(2*nAVG + 1);
for quasar_ind = 1:num_quasars
    tic;
    z_qso = z_qsos(quasar_ind);
    fprintf('processing quasar %i/%i (z_QSO = %0.4f) ...', ...
                              quasar_ind, num_quasars, z_qso);
    
    this_wavelengths    =    all_wavelengths{quasar_ind};
    this_wavelengths    =              this_wavelengths';
    this_flux           =           all_flux{quasar_ind}; 
    this_flux           =                     this_flux';
    this_noise_variance = all_noise_variance{quasar_ind};
    this_noise_variance =           this_noise_variance';
    this_pixel_mask     =     all_pixel_mask{quasar_ind};
    this_pixel_mask     =               this_pixel_mask';
    % this_sigma_pixel    =     all_sigma_pixel{quasar_ind};
    % this_sigma_pixel    =               this_sigma_pixel';
    
    % convert to QSO rest frame
    this_rest_wavelengths = emitted_wavelengths(this_wavelengths, z_qso);
    
    unmasked_ind = (this_rest_wavelengths >= min_lambda) & ...
        (this_rest_wavelengths <= max_lambda);% & (this_sigma_pixel>0);
    % keep complete copy of equally spaced wavelengths for absorption
    % computation
    this_unmasked_wavelengths = this_wavelengths(unmasked_ind);
    
    % [mask_ind] remove flux pixels with pixel_mask; pixel_mask is defined
    % in read_spec_dr7.m
    ind = unmasked_ind & (~this_pixel_mask);% & (this_sigma_pixel>0);
    
    this_wavelengths      =      this_wavelengths(ind);
    this_rest_wavelengths = this_rest_wavelengths(ind);
    this_flux             =             this_flux(ind);
    this_noise_variance   =   this_noise_variance(ind);
    % this_sigma_pixel      =      this_sigma_pixel(ind);
    %   this_lya_zs = ...
    %       (this_wavelengths - lya_wavelength) / ...
    %       lya_wavelength;
    
    % c4 existence prior
    less_ind = (prior.z_qsos < (z_qso + prior_z_qso_increase));
    this_num_c4    = nnz(prior.c4_ind(less_ind));
    this_num_quasars = nnz(less_ind);
    this_p_c4 = this_num_c4 / this_num_quasars;
    log_priors_c4(quasar_ind) = ...
        log(                   this_num_c4) - log(this_num_quasars);
    log_priors_no_c4(quasar_ind) = ...
        log(this_num_quasars - this_num_c4) - log(this_num_quasars);
    fprintf('\n');
    fprintf(' ...     p(   CIV | z_QSO)  mvn      : %0.3f\n',     this_p_c4);
    fprintf(' ...     p(no CIV | z_QSO)        : %0.3f\n', 1 - this_p_c4);
    
    % interpolate model onto given wavelengths
    this_mu = mu_interpolator( this_rest_wavelengths);
    this_M  =  M_interpolator({this_rest_wavelengths, 1:k});
    
    %   this_log_omega = log_omega_interpolator(this_rest_wavelengths);
    %   this_omega2 = exp(2 * this_log_omega);
    
    %   this_scaling_factor = 1 - exp(-tau_0 .* (1 + this_lya_zs).^beta) + c_0;
    
    %   this_omega2 = this_omega2 .* thiscaling_factor.^2;
    
    % baseline: probability of no DLA model
    %   disp(size(this_M))
    %   disp(size(this_flux))
    %   disp(size(this_noise_variance))
    log_likelihoods_no_c4(quasar_ind) = ...
        log_mvnpdf_low_rank(this_flux, this_mu, this_M, ...
        this_noise_variance);
    
    log_posteriors_no_c4(quasar_ind) = ...
        log_priors_no_c4(quasar_ind) + log_likelihoods_no_c4(quasar_ind);
    
    fprintf(' ... log p(D | z_QSO, no CIV)     : %0.2f\n', ...
        log_likelihoods_no_c4(quasar_ind));
    
    min_z_c4s(quasar_ind) = min_z_c4(this_wavelengths, z_qso);
    max_z_c4s(quasar_ind) = max_z_c4(z_qso, max_z_cut);
    
    sample_z_c4 = ...
        min_z_c4s(quasar_ind) +  ...
        (max_z_c4s(quasar_ind) - min_z_c4s(quasar_ind)) * offset_z_samples;


    % Temperature samples
    
    % ensure enough pixels are on either side for convolving with
    % instrument profile

   % building a finer wavelength and mask arrays 
   % by adding the mean of ith and ith +1 element
    lenW = length(this_wavelengths);
    % fine_wavelength = zeros(2*lenW-1,1);
    % fine_pixel_mask = zeros(2*lenW-1,1);
    % for i=1:lenW-1
        
    %     fine_wavelengths(2*i) = (this_wavelengths(i) + this_wavelengths(i+1))/2;
    %     fine_wavelengths(2*i-1) = this_wavelengths(i);
    %     fine_pixel_mask(2*i) = (this_pixel_mask(i) + this_pixel_mask(i+1))/2;
    %     fine_pixel_mask(2*i-1) = this_pixel_mask(i);
        

    % end
    % fine_wavelengths(2*lenW-1) = this_wavelengths(end);
    % fine_pixel_mask(2*lenW-1) = this_pixel_mask(end);
    fine_wavelengths = finer(this_wavelengths, nAVG);
    fine_pixel_mask = finer(this_pixel_mask, nAVG);
    fine_pixel_mask(fine_pixel_mask>0)=1;

    this_rest_fine_wavelengths = emitted_wavelengths(fine_wavelengths, z_qso);
    unmasked_ind_fine = (this_rest_fine_wavelengths >= min_lambda) & ...
    (this_rest_fine_wavelengths <= max_lambda);
    this_unmasked_fine_wavelengths = fine_wavelengths(unmasked_ind_fine);

    padded_fine_wavelengths = ...
        [logspace(log10(min(this_unmasked_fine_wavelengths)) - width * pixel_spacing/2, ...
        log10(min(this_unmasked_fine_wavelengths)) - pixel_spacing/2,...
        width)';...
        this_unmasked_fine_wavelengths';...
        logspace(log10(max(this_unmasked_fine_wavelengths)) + pixel_spacing/2,...
        log10(max(this_unmasked_fine_wavelengths)) + width * pixel_spacing/2,...
        width)'...
        ];
    
    % [mask_ind] to retain only unmasked pixels from computed absorption profile
    % this has to be done by using the unmasked_ind which has not yet
    % been applied this_pixel_mask.
    ind = (~this_pixel_mask(unmasked_ind));
    
    % compute probabilities under DLA model for each of the sampled
    % (normalized offset, log(N HI)) pairs
    lenW_unmasked = length(this_unmasked_wavelengths);
    parfor i = 1:num_C4_samples
        
        % Limitting red-shift in the samples
        this_z_c4 = (this_wavelengths / 1549) - 1;
        % absorption corresponding to this sample with one absorption line as a noise model 

        % num_lines=1;
        % absorptionL1 = voigt(padded_wavelengths, sample_z_c4(i), ...
        %     nciv_samples(i), num_lines, sample_sigma_c4(i), this_sigma_pixel);
        
        % absorptionL1 = absorptionL1(ind);
        % c4_muL1     = this_mu     .* absorptionL1;
        % c4_ML1      = this_M      .* absorptionL1;
        % %     dla_omega2 = this_omega2 .* absorption.^2;
        
        % sample_log_likelihoods_c4L1(quasar_ind, i) = ...
        %     log_mvnpdf_low_rank(this_flux, c4_muL1, c4_ML1, ...
        %     this_noise_variance);

        % absorption corresponding to this sample with two absorption lines as a doublet 
        % compute fine absorption 
        num_lines=2;
        absorptionL2_fine = voigt0(padded_fine_wavelengths, sample_z_c4(i), ...
        nciv_samples(i),num_lines , sample_sigma_c4(i));
        
        % average fine absorption and shrink it to the size of original array
        % as large as the unmasked_wavelengths
        
        absorptionL2 = zeros(lenW_unmasked, 1);
        absorptionL2(1,1) = absorptionL2_fine(1);
        absorptionL2(end,1) = absorptionL2_fine(end);
        for ii=2:lenW_unmasked-1
            s=0;
            for j=0:2*nAVG
                s=s+absorptionL2_fine((nAVG+1)*ii-j);
            end
                absorptionL2(ii,1) = s*fAVG;
        end

        absorptionL2 = absorptionL2(ind);
        c4_muL2     = this_mu     .* absorptionL2;
        c4_ML2      = this_M      .* absorptionL2;
        %     dla_omega2 = this_omega2 .* absorption.^2;
        sample_log_likelihoods_c4L2(quasar_ind, i) = ...
          log_mvnpdf_low_rank(this_flux, c4_muL2, c4_ML2, ...
          this_noise_variance);
    %       figure()
    %       plot(this_unmasked_wavelengths, absorptionL2, '-o', 'Color', 'r')
    %       hold on
    %       plot(fine_wavelengths, absorptionL2_fine, '-o', 'Color', 'k')
    end
       
    % compute sample probabilities and log likelihood of DLA model in
    % numerically safe manner for one line
    % % max_log_likelihoodL1 = max(sample_log_likelihoods_c4L1(quasar_ind, :));
    % % sample_probabilitiesL1 = ...
    % %     exp(sample_log_likelihoods_c4L1(quasar_ind, :) - ...
    % %     max_log_likelihoodL1);
    % % log_likelihoods_c4L1(quasar_ind) = ...
    % %     max_log_likelihoodL1 + log(mean(sample_probabilitiesL1));
    
    % % log_posteriors_c4L1(quasar_ind) = ...
    % %     log_priors_c4(quasar_ind) + log_likelihoods_c4L1(quasar_ind);
    
    % fprintf(' ... log p(D | z_QSO,    L1)     : %0.2f\n', ...
    %     log_likelihoods_c4L1(quasar_ind));
    % fprintf(' ... log p(L1 | D, z_QSO)        : %0.2f\n', ...
    %     log_posteriors_c4L1(quasar_ind));

       
    % compute sample probabilities and log likelihood of DLA model in
    % numerically safe manner for  doublet 
    max_log_likelihoodL2 = max(sample_log_likelihoods_c4L2(quasar_ind, :));
    sample_probabilitiesL2 = ...
        exp(sample_log_likelihoods_c4L2(quasar_ind, :) - ...
        max_log_likelihoodL2);
    log_likelihoods_c4L2(quasar_ind) = ...
        max_log_likelihoodL2 + log(mean(sample_probabilitiesL2));
    
    log_posteriors_c4L2(quasar_ind) = ...
        log_priors_c4(quasar_ind) + log_likelihoods_c4L2(quasar_ind);
    
    fprintf(' ... log p(D | z_QSO,    CIV)     : %0.2f\n', ...
        log_likelihoods_c4L2(quasar_ind));
    fprintf(' ... log p(CIV | D, z_QSO)        : %0.2f\n', ...
        log_posteriors_c4L2(quasar_ind));
    %  fprintf(' ... Num_CIV                      : %d\n ', ...
    %    Full_catalog.all_Num_c4_sys(quasar_ind))
    % fprintf('... FilterFlag                    : %d\n ', filter)
    
    % [~, maxindL1] = nanmax(sample_log_likelihoods_c4L1(quasar_ind, :));
    % map_z_c4L1(quasar_ind)    = sample_z_c4(maxindL1);        
    % map_N_c4L1(quasar_ind)  = log_nciv_samples(maxindL1);
    % map_sigma_c4L1(quasar_ind)  = sample_sigma_c4(maxindL1);
    % fprintf('L1\nmap(N): %.2f, map(z_c4): %.2f, map(b/1e5): %.2f\n',map_N_c4L1(quasar_ind),...
    %  map_z_c4L1(quasar_ind), map_sigma_c4L1(quasar_ind)/1e5);

    [~, maxindL2] = nanmax(sample_log_likelihoods_c4L2(quasar_ind, :));
    map_z_c4L2(quasar_ind)    = sample_z_c4(maxindL2);        
    map_N_c4L2(quasar_ind)  = log_nciv_samples(maxindL2);
    map_sigma_c4L2(quasar_ind)  = sample_sigma_c4(maxindL2);
    fprintf('L2\nmap(N): %.2f, map(z_c4): %.2f, map(b/1e5): %.2f\n',map_N_c4L2(quasar_ind),...
    map_z_c4L2(quasar_ind), map_sigma_c4L2(quasar_ind)/1e5);
    fprintf(' took %0.3fs.\n', toc);
    
end

% compute model posteriors in numerically safe manner
max_log_posteriors = ...
    max([log_posteriors_no_c4,  log_posteriors_c4L2], [], 2);

model_posteriors = ...
    exp([log_posteriors_no_c4,  log_posteriors_c4L2] - max_log_posteriors);

model_posteriors = model_posteriors ./ sum(model_posteriors, 2);

p_no_c4 = model_posteriors(:, 1);

p_c4    = 1 - p_no_c4;

% save results
variables_to_save = {'release', 'training_set_name', ...
     'prior_ind', 'release', ...
    'test_ind', 'prior_z_qso_increase', ...
    'max_z_cut', 'min_z_c4s', 'max_z_c4s', ...
    'log_priors_no_c4', 'log_priors_c4', ...
    'log_likelihoods_no_c4',  ...
    'sample_log_likelihoods_c4L2', 'log_likelihoods_c4L2'...
    'log_posteriors_no_c4', 'log_posteriors_c4L1', 'log_posteriors_c4L2',...
    'model_posteriors', 'p_no_c4', ...
    'map_z_c4L2', 'map_N_c4L2', 'map_sigma_c4L2', 'p_c4'};

filename = sprintf('%s/processed_qsos_R%s', ...
    processed_directory(release), ...
    testing_set_name);

save(filename, variables_to_save{:}, '-v7.3');