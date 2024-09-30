%     Copyright (C) 2019  Emanuele Ortu, Eleonora Sulas, Danilo Pani.
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.

% To contact the authors:
% sulaseleonora1992@gmail.com, danilo.pani@unica.it,
% Address:  Department of Electrical and Electronic Engineering, 
% University of Cagliari,Via Marengo 3 - 09123 Cagliari, Italy

% Cagliari, 15 Jan 2020

% Please cite the following article when using this code:
% TBA

% This matlab code extract the upper and the lower envelope of the Doppler 
% image.
% input:
%   filename: the filename as a string and bmp format
% output:
%   x_up: upper envelope
%   x_down: lower envelope

function [x_up,x_down]=envelope_extraction(filename)

    if nargin==0
        [file,path] = uigetfile({'*.bmp'});
        filename=fullfile(path,file);
    end
    % loading
    x = imread(filename);
    x = rgb2gray(x);
    % delete x axis line(from white to black)
    ascissa = find(sum(x')==255*length(x));
    x(ascissa,:)= 0;    
    % median otsu 2D threshold
    T = median_otsu2d(x,5);   
    
    BW = imbinarize(x,T/255);
    BW = bwareaopen(BW,70,4);
 
    x1 = BW;
    x_down = [];
    x_up = [];
    up = x1(1:ascissa-1,:);
    down = x1(ascissa+1:end,:);
    
    for j = 1:length(x1)        
        
        ind = find(down(:,j)==1);
        
        if isempty(ind) == 1
            x_down(j) = 0;
        else
            x_down(j)= ind(end); 
        end
        
        ind = find(up(:,j)==1);
        if isempty(ind) == 1
            x_up(j) = ascissa;
        else
            x_up(j)= ind(1); 
        end
    end
    
    x_down = - x_down;
    x_up = ascissa - (x_up);   
          
    
end

function th = median_otsu2d(img1,k)

dim = size(img1);
if length(dim)==3
    img1 = rgb2gray(img1);
end

k_2=fix(k/2);
[m,n]=size(img1);

f=img1(1+k_2:m-k_2,1+k_2:n-k_2);
g=zeros(m-k_2,n-k_2);
f=double(f);
c=zeros(256,256);

for i=1+k_2:m-k_2
    for j=1+k_2:n-k_2
        g(i-k_2,j-k_2)=median(median((img1(i-k_2:i+k_2,j-k_2:j+k_2))));  
        c(f(i-k_2,j-k_2)+1,g(i-k_2,j-k_2)+1)=c(f(i-k_2,j-k_2)+1,g(i-k_2,j-k_2)+1)+1;
    end
end
[m1,n1]=size(f);
p = c./(m1*n1);
figure
surf(p)

for s=0:255
    for t=0:255        
       w0=sum(sum(p(1:s,1:t)));
       w1=sum(sum(p(s+1:end,t+1:end)));
       u0=[sum(sum(p(1:s,1:t)'*(1:s)'))/w0, sum(sum(p(1:s,1:t)*(1:t)'))/w0 ];
       u1=[sum(sum(p(s+1:end,t+1:end)'*(s+1:256)'))/w1, sum(sum(p(s+1:end,t+1:end)*(t+1:256)'))/w1];
       uT=[sum(sum(p(1:end,1:end)'*(1:256)')), sum(sum(p(1:end,1:end)*(1:256)'))];
       sigma=(w0*((u0-uT)*(u0-uT)'))+(w1*((u1-uT)*(u1-uT)'));
       trace(s+1,t+1)=(w0*(((u0(1)-uT(1))^2)+(u0(2)-uT(2))^2))+(w1*(((u1(1)-uT(1))^2)+(u1(2)-uT(2))^2));      
    end
end
[mas,ind] = max(trace);
[mas2, ind2] = max(mas);
th = ind2;
end

%Code Start
%run envelope extraction
[x_up, x_down] = envelope_extraction('DataSetForAssignment.bmp');

%Calculate heart rate using autocorrelation
function heart_rate = calculate_heart_rate(envelope, frequency)

    [acf, lags] = xcorr(envelope, 'coeff');
    
    [~, maxLag] = max(acf(lags > 0)); 
    peakLag = lags(maxLag + find(lags == 0));
    
    period = peakLag / frequency;
    heart_rate = 60 / period;
end

function heart_rate_series = sliding_window(envelope, frequency, window_size, step_size)
    samples = round(window_size * frequency);
    step_samples = round(step_size * frequency);
    num = floor((length(envelope) - samples) / step_samples) + 1;
    
    heart_rate_series = zeros(1, num);
    
    for i = 1:num
        start_idx = (i-1) * step_samples + 1;
        end_idx = start_idx + samples - 1;
        window_data = envelope(start_idx:end_idx);
        
        heart_rate_series(i) = calculate_heart_rate(window_data, frequency);
    end
end

frequency = 1000; 
window = 3.75;  
increment = 0.75; 

%calculate heart rate series and time series
heart_rate_series = sliding_window(x_up, frequency, window, increment);
time_series = (0:(length(heart_rate_series)-1)) * 0.75;  


% Plot heart rate over time
figure;
plot(time_series, heart_rate_series, 'LineWidth', 3);
xlabel('Time (s)');
ylabel('Heart Rate (BPM)');
title('Heart Rate vs Time');

%parts of code understood and generate with internet search and GPT.