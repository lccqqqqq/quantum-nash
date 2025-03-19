% creating 3D plots for the SU(2) matrices
data = zeros(4,1500);

data(1,:) = reshape(real(GameStrategy{2}(1,1,:)),1,1500); % X input, Pauli I
data(2,:) = reshape(imag(GameStrategy{2}(1,1,:)),1,1500); % Y input, Pauli iZ
data(3,:) = reshape(real(GameStrategy{2}(1,1,:)),1,1500); % Z input, Pauli X
data(4,:) = reshape(imag(GameStrategy{2}(1,1,:)),1,1500); % U input, Pauli -Y

% Stereographic projection to R3

data_sp = zeros(3,1500); % for stereographically projected coordinates
data_sp(1,:) = 2*data(2,:)./(1+data(1,:));
data_sp(2,:) = 2*data(3,:)./(1+data(1,:));
data_sp(3,:) = 2*data(4,:)./(1+data(1,:)); 

% plot the projected path in R3
plot3(data_sp(1,1:1500),data_sp(2,1:1500),data_sp(3,1:1500),'r.')