function [ ET_serie ET ] = ET_Thornthwaite( T_serie , tdates , lat  )
% [ ET_serie ] = ET_Thornthwaite( T_serie , tdates , lat  )
% Compute evapotranspiration (ET) using the ET_Thornthwaite method.
% References:
%   - Thornthwaite CW. 1948. An approach toward a rational classification 
%     of climate. Geographical Review 38: 55–94.
%   - Xu C.Y. and Singh V.P. 2001.  Evaluation and generalization of 
%     temperature-based methods for calculating evaporation. 
%     Hydrological Processes 15, 305–319.
% 
% Input:
%   - T_serie    = time serie of temperature   [mm]   - vector [ h , 1 ]
%   - tdates     = time [dd,mm,yyyy]                  - matrix [ h , 3 ]
%   - lat        = North latitude              [°]    - scalar
%                  acceptable values = [0:10:30 , 35:5:50] 
% Output:
%   - ET_serie   = time serie of montly mean 
%                  evapotranspiration          [mm]   - vector [ h , 1 ]
% 
% NOTE: first data in T_serie must refer to the 1st of January and last
% data must refer to the 31th of December
% 
% Last update: Daniela Anghileri 02/08/2014

% input control
if ~any([0:10:30 , 35:5:50] == lat)
    error( 'Error in ET_Thornthwaite: latitude must be one of the following values: 0    10    20    30    35    40    45    50') ;
end

% all days but 29s of february
idx = ( tdates( : , 1 ) .* tdates( : , 2) ) ~= 58 ; % because 58= 29 x II (February)
T   = T_serie( idx ) ;   % time serie without leap day
N   = length(T)/365 ;       % number of years

% compute montly mean
A = reshape( T , 365 , N )'                                  ;
%                     J  F  M  A  M  J  J  A  S  O  N  D
day_per_month     = [ 31 28 31 30 31 30 31 31 30 31 30 31 ] ;
day_per_month_cum = [ 0 cumsum(day_per_month) ]             ;
Tmm = NaN( N , 12 ) ;
    for j = 1:12
        mm  = mean( A( : , day_per_month_cum(j)+1 : day_per_month_cum(j+1) ) , 2 );
        Tmm(:,j) = mm  ;
    end

% compute Thornthwaite formula
    i = (Tmm/5).^1.51 ;
    i( Tmm<=0 ) = 0   ;
    I = sum( i , 2 )  ;
    a = 67.5 *10^-8 * I.^3 - 77.1 *10^-6 * I.^2 + 0.0179 * I + 0.492         ;
% correction factor depending on latitute
%   latN  J    F    M    A    M    J    J    A    S    O    N    D
    BS = [ 0 1.04 0.94 1.04 1.01 1.04 1.01 1.04 1.04 1.01 1.04 1.01 1.04 
          10 1.00 0.91 1.03 1.03 1.08 1.06 1.08 1.07 1.02 1.02 0.98 0.99
          20 0.95 0.90 1.03 1.05 1.13 1.11 1.14 1.11 1.02 1.00 0.93 0.94
          30 0.90 0.87 1.03 1.08 1.18 1.17 1.20 1.14 1.03 0.98 0.89 0.88
          35 0.87 0.85 1.03 1.09 1.21 1.21 1.23 1.16 1.03 0.97 0.86 0.85
          40 0.84 0.83 1.03 1.11 1.24 1.25 1.27 1.18 1.04 0.96 0.83 0.81
          45 0.80 0.81 1.02 1.13 1.28 1.29 1.31 1.21 1.04 0.94 0.79 0.75 
          50 0.74 0.78 1.02 1.15 1.33 1.36 1.37 1.25 1.06 0.92 0.76 0.70 ] ;
    b = BS( BS(:,1)==lat , 2:end )       ; % [1,12]
    ET = 16 * repmat(b,N,1) .* ( 10*Tmm ./ repmat(I,1,12) ) .^ repmat(a,1,12) ;
    ET( Tmm<=0 ) = 0 ;
    
    %ET_month  = reshape( ET' , [] , 1 )   ;
    %ET_monthly_mean=mean(ET,1);
    
% % VARIANTE
% % "... Volendo estenderne l’uso alla stima dell’evapotraspirazione 
% % potenziale dei singoli mesi di una serie storica, occorre ovviamente 
% % assumere I costante per una data località, uguale al valore 
% % corrispondente alle medie delle temperature mensili calcolate su un 
% % periodo abbastanza lungo..." (fonte: dispense di idrologia trovate in
% % internet)
% Tmm_mean_year = mean( Tmm )     ;    % anno medio: [1,12]
% i = (Tmm_mean_year/5).^1.51     ;
% i( Tmm_mean_year<=0 ) = 0       ;
% I = sum( i ) ;
% % correction factor depending on latitute
% %   latN  J    F    M    A    M    J    J    A    S    O    N    D
% BS = [ 0 1.04 0.94 1.04 1.01 1.04 1.01 1.04 1.04 1.01 1.04 1.01 1.04 
%       10 1.00 0.91 1.03 1.03 1.08 1.06 1.08 1.07 1.02 1.02 0.98 0.99
%       20 0.95 0.90 1.03 1.05 1.13 1.11 1.14 1.11 1.02 1.00 0.93 0.94
%       30 0.90 0.87 1.03 1.08 1.18 1.17 1.20 1.14 1.03 0.98 0.89 0.88
%       35 0.87 0.85 1.03 1.09 1.21 1.21 1.23 1.16 1.03 0.97 0.86 0.85
%       40 0.84 0.83 1.03 1.11 1.24 1.25 1.27 1.18 1.04 0.96 0.83 0.81
%       45 0.80 0.81 1.02 1.13 1.28 1.29 1.31 1.21 1.04 0.94 0.79 0.75 
%       50 0.74 0.78 1.02 1.15 1.33 1.36 1.37 1.25 1.06 0.92 0.76 0.70 ] ;
% 
% b = BS( BS(:,1)==lat , 2:end )       ; % [1,12]
% a = 67.5 *10^-8 * I^3 - 77.1 *10^-6 * I^2 + 0.0179 * I + 0.492         ;
% ET = 16.2 * repmat( b , N , 1) .* ( 10/I * Tmm ) .^ a ; % [N,12]  
% ET( Tmm<=0 ) = 0 ;
  

% compute ET for each day of the serie
    C = NaN( N , 365 ) ;
        for j =1:12
            %c = repmat( ET(:,j) , 1 , day_per_month(j) ) ;     % [N,day_per_month]
            c = repmat( ET(:,j)/day_per_month(j) , 1 , day_per_month(j) ) ;     % [N,day_per_month]
            C( : , day_per_month_cum(j)+1 : day_per_month_cum(j+1) ) = c ;
        end
    ET_serie       = NaN( size(T_serie) )     ;
    ET_serie(idx)  = reshape( C' , [] , 1 )   ;
    ET_serie(~idx) = ET_serie( find(~idx)-1 ) ;
