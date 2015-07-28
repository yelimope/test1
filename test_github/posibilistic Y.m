%% Possibilistic Consistency Analysis of datasets
 
clear all  % load and prepare
 
%opening and load Matriz N

[DATA]     =xlsread('MODEL_matrix_2012.xls','', '', 'basic'); % Opening model
irrev      = not(DATA(3,1:end));
N          = DATA(4:40,1:end);
N(isnan(N))=0;
clear DATA; 

%opening and load the experimental data

[DATA,TXT]=xlsread('DATA_measurements_actualizados _2.xls','', '', '');
DATA = DATA(1:end,:);

%% Perform PMFA computations

DATOS  = [];
for i=1:size(DATA,1)
    
	%  measures  OUR   			GLU     	CER     	Et      Gly   		  Cit    	 Pyr            MET         Bio   		     
    ED(i).ind = [ 40             41          42         43       44            45         46             47         48];
    ED(i).wm = [DATA(i,14)   DATA(i,8)   DATA(i,15)   DATA(i,13)   DATA(i,9)  DATA(i,12)   DATA(i,11)  DATA(i,10)  DATA(i,7)	]';
	ED(i).protein = DATA(i,18);
    ED(i).biomass = DATA(i,7);
    ED(i).indice = DATA(i,1);
    
	% Remove measurements NAN
	ED(i).ind(isnan(ED(i).wm)) = [];  ED(i).protein(isnan(ED(i).wm))=[]; ED(i).wm(isnan(ED(i).wm)) = [];
    
   % Constraint-based model
    [CBM] = POS_define_MOC(N, irrev);
    [CBM] = POS_define_MEC(CBM, ED(i), 0.05, 0.2, 0.001);   % this is the criteria for uncertainty (see Tortajada 2012)
    
	%% most possible solution    
	solvesdp(CBM.CB,CBM.J);
    poss = exp(-double(CBM.J));
    
    Fluxmp  = double(CBM.v(45));
    
    %% Use interval
    Fluxp1  = POS_interval(CBM.CB, CBM.J, CBM.v(45), 0.99, 'cond');
    DATOS  = [DATOS; ED(i).indice Fluxmp ED(i).biomass ED(i).protein Fluxp1];
    
    i
end

%%
DATA_MFA = DATOS;
save DATA_MFA DATA_MFA