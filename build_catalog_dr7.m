% Build catalogs usable for spectra from dr16

Seyffert_MgII_detected = fitsread(...
'data/MGII_catalogs/Seyffert_MGII_cat/distfiles/sdss_mgii_seyffertetal13.fit',...
'binarytable');

mgII_QSO_ID                   = Seyffert_MgII_detected{1};
z_qso_system                  = Seyffert_MgII_detected{10};
Z_abs_ORG                     = Seyffert_MgII_detected{17};
EW                            = Seyffert_MgII_detected{22};
SigmaEW                       = Seyffert_MgII_detected{23};
flagEW                        = Seyffert_MgII_detected{24};
Column_density_MGII_ORG       = Seyffert_MgII_detected{27};
ERR_Column_density_MGII_ORG   = Seyffert_MgII_detected{28};
Column_density_MGII_FLG       = Seyffert_MgII_detected{29};
dummy                         = Seyffert_MgII_detected{30};
RATING                        = dummy(:,1);
% % filtering out those column densities with not good measurements
EW1                           = EW(:,1);
EW2                           = EW(:,2);
[nSys,dd]=size(mgII_QSO_ID);
NMGII=zeros(nSys,1);
Z_mgii=zeros(nSys,1);
for i=1:nSys
    NMGII(i) = Column_density_MGII_ORG(i,1)/ERR_Column_density_MGII_ORG(i,1)^2 + Column_density_MGII_ORG(i,2)/ERR_Column_density_MGII_ORG(i,2)^2;
    NMGII(i)=NMGII(i)/(1/ERR_Column_density_MGII_ORG(i,1)^2+1/ERR_Column_density_MGII_ORG(i,2)^2);
    Z_mgii(i) = min(Z_abs_ORG(i,1),Z_abs_ORG(i,2));
end

save('data/MGII_catalogs/Seyffert_MGII_cat/processed/MGII-cat.mat','mgII_QSO_ID','Z_mgii','NMGII');

% There are some NAN valued c4_NCIV
% extract basic QSO information from Cookse_all_QSO catalog 
quasar_catalog = ...
fitsread('data/dr7/distfiles/dr7qso_MgII_noBAL.fit', 'binarytable');
all_plate_dr7             = quasar_catalog{48};
all_mjd_dr7             = quasar_catalog{47};
all_fiber_dr7             = quasar_catalog{49};
all_RA                = quasar_catalog{2};
all_DEC               = quasar_catalog{3};
all_zqso                = quasar_catalog{4};
num_quasars             = numel(all_zqso);
all_QSO_ID=cell(num_quasars,1);

all_z_mgii = zeros(num_quasars, 20)-1;
all_N_mgii = zeros(num_quasars, 20)-1;
all_RATING = zeros(num_quasars, 20)-1;
all_EW1 = zeros(num_quasars, 20)-1;
all_EW2 = zeros(num_quasars, 20)-1;
for i=1:num_quasars
    all_QSO_ID{i}=sprintf('%05i-%04i-%03i', (all_mjd_dr7(i)), ...
    (all_plate_dr7(i)), (all_fiber_dr7(i)));
    ThisSystems = ismember(mgII_QSO_ID, all_QSO_ID{i});
    thisZ_mgIIs = Z_mgii(ThisSystems);
    thisN_mgIIs = NMGII(ThisSystems);
    this_RATING = RATING(ThisSystems);
    this_EW1 = EW1(ThisSystems);
    this_EW2 = EW1(ThisSystems);
    nSys = nnz(ThisSystems);
    
    for j=1:nSys
        all_z_mgii(i,j) = thisZ_mgIIs(j);
        all_N_mgii(i,j) = thisN_mgIIs(j);
        all_RATING(i,j) = this_RATING(j);
        all_EW1(i, j) = this_EW1(j);
        all_EW2(i, j) = this_EW2(j);

    end
end



% save catalog 
% need to fix directory here?
release = 'dr7';
variables_to_save = {'all_plate_dr7', 'all_mjd_dr7', 'all_fiber_dr7', ...
 'all_QSO_ID', 'all_RA', 'all_DEC', 'all_zqso', 'all_EW1', 'all_EW2',...
 'all_N_mgii','all_z_mgii', 'all_RATING', 'mgII_QSO_ID'};
save(sprintf('%s/catalog', processed_directory(release)), ...
    variables_to_save{:}, '-v7.3');
