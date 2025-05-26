clc;
clear;
close all;


%% Učitavanje slike (Loading the image)
% Učitavanje originalne slike i njeno normalizovanje u opseg [0, 1]
imagePath = 'C:\Users\Bato\Desktop\RADOVI 2025\Pigmentacija\img.png';
originalImage = imread(imagePath);
originalImage = im2double(originalImage);

% Provera da li slika ima 3 kanala (R, G, B)
[rows, cols, channels] = size(originalImage);
if channels ~= 3
    error('Slika mora imati 3 kanala (R, G, B) / Image must have 3 channels (R, G, B)');
end

%% Ekstrakcija R, G, B kanala (Extracting R, G, B channels)
% Izdvajanje crvenog (R), zelenog (G) i plavog (B) kanala
R = originalImage(:,:,1);
G = originalImage(:,:,2);
B = originalImage(:,:,3);

%% Pojačavanje nijansi uz Lorencovu krivu (Enhancing hues using Lorenz curve)
% Kreiranje Lorencove krive za nijansiranje boja
x = linspace(0, 1, max(rows, cols)); % Diskretno polje
lorenzCurve = 1 ./ (1 + (x - 0.5).^2); % Lorencova kriva

% Faktor intenziteta za prilagođavanje nijansi
scaleFactor = 0.05; % Minimalne promene boja

% Pojačavanje prelaza između boja (Enhancing color transitions)
RG_mix = (1 - scaleFactor) * R + scaleFactor * (R .* (1 - lorenzCurve') + G .* lorenzCurve'); % R -> G
GB_mix = (1 - scaleFactor) * G + scaleFactor * (G .* (1 - lorenzCurve') + B .* lorenzCurve'); % G -> B
BR_mix = (1 - scaleFactor) * B + scaleFactor * (B .* (1 - lorenzCurve') + R .* lorenzCurve'); % B -> R

% Kreiranje novih kanala sa pojačanim nijansama
R_enhanced = RG_mix .* (1 - GB_mix); % Dominacija crvene (R dominance)
G_enhanced = GB_mix .* (1 - BR_mix); % Dominacija zelene (G dominance)
B_enhanced = BR_mix .* (1 - RG_mix); % Dominacija plave (B dominance)

% Kombinovanje pojačanih kanala u jednu sliku
enhancedImage = cat(3, R_enhanced, G_enhanced, B_enhanced);

%% Interpolacija sa originalnom slikom (Interpolation with the original image)
% Maksimalna dozvoljena razlika u nijansama
maxDifference = 0.3; % 30% maksimalne razlike
finalImage = (1 - maxDifference) * originalImage + maxDifference * enhancedImage;

%% Posvetljavanje i dodatno potamnjivanje tamnih tonova (Brightening and darkening)
% Faktor posvetljavanja
brightnessFactor = 1.4; % Povećanje svetline (svetliji tonovi)
brightenedImage = finalImage * brightnessFactor; % Povećanje osvetljenja

% Gama korekcija za nežnije tonove
brightenedImage = brightenedImage .^ 0.8;

% Potamnjivanje tamnih tonova za 7%
darkenFactor = 0.9; % Smanjenje za 7%
darkPixelThreshold = 0.5; % Prag za tamne piksele
darkenedImage = brightenedImage;
darkenedImage(brightenedImage < darkPixelThreshold) = ...
    darkenedImage(brightenedImage < darkPixelThreshold) * darkenFactor;

% Normalizacija osvetljenja u opseg [0, 1]
darkenedImage = min(darkenedImage, 1); % Ograničavanje na maksimalnu vrednost 1

%% Generisanje čistih boja za R, G i B (Generating pure R, G, and B channels)
R_pure = cat(3, R, zeros(size(G)), zeros(size(B))); % Čista crvena (Pure red)
G_pure = cat(3, zeros(size(R)), G, zeros(size(B))); % Čista zelena (Pure green)
B_pure = cat(3, zeros(size(R)), zeros(size(G)), B); % Čista plava (Pure blue)

%% Prikaz rezultata (Displaying results)

% Prikaz originalne slike (Displaying original image)
figure; imshow(originalImage); title('Originalna slika / Original Image');

% Prikaz posvetljene i potamnjene slike (Displaying brightened and darkened image)
figure; imshow(darkenedImage); title('Posvetljena i potamnjena slika / Brightened and Darkened Image');
imwrite(darkenedImage,"img-LPM.png")

% Prikaz čistih boja (Displaying pure colors)
figure; imshow(R_pure); title('Čista crvena (R) / Pure Red (R)');
figure; imshow(G_pure); title('Čista zelena (G) / Pure Green (G)');
figure; imshow(B_pure); title('Čista plava (B) / Pure Blue (B)');

%% Prikaz 3D grafa intenziteta (3D intensity surface plot)
[xGrid, yGrid] = meshgrid(1:cols, 1:rows);
intensity = sum(darkenedImage, 3); % Intenzitet boje (Color intensity)
figure;
surf(xGrid, yGrid, intensity, 'EdgeColor', 'none');
colormap jet; colorbar; title('3D grafikon intenziteta / 3D Intensity Surface Plot');

%% Histogram razlika u skali sive (Histogram of differences in grayscale)
difference = abs(originalImage - darkenedImage); % Razlika između slika
figure;
histogram(difference(:), 50, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'k'); % Histogram sa sivim tonovima
colormap gray; % Skala sive boje
title('Histogram razlika / Histogram of Differences');

%% Izračunavanje SSIM (Calculating SSIM)
[ssimVal, ssimMap] = ssim(darkenedImage, originalImage); % Izračunavanje SSIM
disp(['SSIM vrednost: ', num2str(ssimVal)]); % Prikaz SSIM vrednosti u komandnom prozoru
figure; imshow(ssimMap, []); title('SSIM mapa / SSIM Map');

%% Poređenje tradicionalnog i novog metoda (Comparison of traditional and new method)

% Razlika u intenzitetu (Intensity difference)
intensityDifference = abs(sum(originalImage, 3) - sum(darkenedImage, 3));
figure;
imagesc(intensityDifference); 
colormap jet; colorbar;
title('Razlika u intenzitetu / Intensity Difference');

% Histogram razlika između tradicionalnog i novog metoda
figure;
histogram(intensityDifference(:), 50, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'k');
colormap gray;
title('Histogram razlika između metoda / Histogram of Differences Between Methods');


% Prelaz R -> G (Transition R → G)
RG_transition = zeros(rows, cols, 3);
RG_transition(:,:,1) = repmat(linspace(1, 0, cols), rows, 1); % Crvena komponenta (Red component)
RG_transition(:,:,2) = repmat(linspace(0, 1, cols), rows, 1); % Zelena komponenta (Green component)
figure; imshow(RG_transition);
title('Prelaz nijansi: R → G / Transition R → G');

% Prelaz R -> B (Transition R → B)
RB_transition = zeros(rows, cols, 3);
RB_transition(:,:,1) = repmat(linspace(1, 0, cols), rows, 1); % Crvena komponenta (Red component)
RB_transition(:,:,3) = repmat(linspace(0, 1, cols), rows, 1); % Plava komponenta (Blue component)
figure; imshow(RB_transition);
title('Prelaz nijansi: R → B / Transition R → B');

% Prelaz G -> B (Transition G → B)
GB_transition = zeros(rows, cols, 3);
GB_transition(:,:,2) = repmat(linspace(1, 0, cols), rows, 1); % Zelena komponenta (Green component)
GB_transition(:,:,3) = repmat(linspace(0, 1, cols), rows, 1); % Plava komponenta (Blue component)
figure; imshow(GB_transition);
title('Prelaz nijansi: G → B / Transition G → B');














%% Generisanje prelaza između kanala uz Lorencovu krivu
% Replikacija Lorencove krive na dimenzije slike (rows x cols)
x = linspace(0, 1, cols); % Diskretno polje duž širine slike
lorenzCurve = 1 ./ (1 + (x - 0.5).^2); % Lorencova kriva
lorenzCurve2D = repmat(lorenzCurve, rows, 1); % Replikacija po redovima

% R -> G prelaz
RG_transition_dif = zeros(rows, cols, 3); % Inicijalizacija prelaza
RG_transition_dif(:,:,1) = R .* (1 - lorenzCurve2D); % Umanjenje crvene komponente
RG_transition_dif(:,:,2) = G .* lorenzCurve2D;       % Povećanje zelene komponente

% G -> B prelaz
GB_transition_dif = zeros(rows, cols, 3); % Inicijalizacija prelaza
GB_transition_dif(:,:,2) = G .* (1 - lorenzCurve2D); % Umanjenje zelene komponente
GB_transition_dif(:,:,3) = B .* lorenzCurve2D;       % Povećanje plave komponente

% B -> R prelaz
BR_transition_dif = zeros(rows, cols, 3); % Inicijalizacija prelaza
BR_transition_dif(:,:,3) = B .* (1 - lorenzCurve2D); % Umanjenje plave komponente
BR_transition_dif(:,:,1) = R .* lorenzCurve2D;       % Povećanje crvene komponente


%% Prikaz 3D surface grafike za prelaze
% R -> G prelaz u 3D
figure(110);
[xGrid, yGrid] = meshgrid(1:cols, 1:rows); % Kreiranje mreže
intensity_RG = RG_transition_dif(:,:,1) + RG_transition_dif(:,:,2); % Zbir boja u prelazu
surf(xGrid, yGrid, intensity_RG, 'EdgeColor', 'none');
colormap jet; colorbar;
title('110. 3D surface prelaza R → G / 3D Surface for R → G Transition');
xlabel('x koordinata / x-coordinate');
ylabel('y koordinata / y-coordinate');
zlabel('Intenzitet prelaza / Transition Intensity');

% G -> B prelaz u 3D
figure(210);
intensity_GB = GB_transition_dif(:,:,2) + GB_transition_dif(:,:,3); % Zbir boja u prelazu
surf(xGrid, yGrid, intensity_GB, 'EdgeColor', 'none');
colormap jet; colorbar;
title('210. 3D surface prelaza G → B / 3D Surface for G → B Transition');
xlabel('x koordinata / x-coordinate');
ylabel('y koordinata / y-coordinate');
zlabel('Intenzitet prelaza / Transition Intensity');

% B -> R prelaz u 3D
figure(310);
intensity_BR = BR_transition_dif(:,:,3) + BR_transition_dif(:,:,1); % Zbir boja u prelazu
surf(xGrid, yGrid, intensity_BR, 'EdgeColor', 'none');
colormap jet; colorbar;
title('310. 3D surface prelaza B → R / 3D Surface for B → R Transition');
xlabel('x koordinata / x-coordinate');
ylabel('y koordinata / y-coordinate');
zlabel('Intenzitet prelaza / Transition Intensity');

%% Prelaz nijansi pomoću naše metode
% Naša metoda koristi diferencijalne transformacije za prelaz
% Primer implementacije za R -> G
alpha = 0.3; % Parametar kontrole prelaza (0-1)

RG_custom = zeros(rows, cols, 3);
RG_custom(:,:,1) = R .* exp(-alpha * lorenzCurve2D); % Umanjenje crvene uz eksponencijalni faktor
RG_custom(:,:,2) = G .* (1 - exp(-alpha * lorenzCurve2D)); % Povećanje zelene

% Prikaz prelaza R -> G uz našu metodu
figure(120);
imshow(RG_custom);
title('120. Prelaz R → G uz našu metodu / R → G Transition with Our Method');

% Implementacija za G -> B
GB_custom = zeros(rows, cols, 3);
GB_custom(:,:,2) = G .* exp(-alpha * lorenzCurve2D); % Umanjenje zelene
GB_custom(:,:,3) = B .* (1 - exp(-alpha * lorenzCurve2D)); % Povećanje plave

% Prikaz prelaza G -> B uz našu metodu
figure(220);
imshow(GB_custom);
title('220. Prelaz G → B uz našu metodu / G → B Transition with Our Method');

% Implementacija za B -> R
BR_custom = zeros(rows, cols, 3);
BR_custom(:,:,3) = B .* exp(-alpha * lorenzCurve2D); % Umanjenje plave
BR_custom(:,:,1) = R .* (1 - exp(-alpha * lorenzCurve2D)); % Povećanje crvene

% Prikaz prelaza B -> R uz našu metodu
figure(320);
imshow(BR_custom);
title('320. Prelaz B → R uz našu metodu / B → R Transition with Our Method');








%% Definicija dimenzija i Lorencove krive
cols = 500; % Broj tačaka u prelazu
x = linspace(0, 1, cols); % Normalizovane pozicije u prelazu
lorenzCurve = 1 ./ (1 + (x - 0.5).^2); % Lorencova kriva

%% Parametri za diferencijalne jednačine
alpha = 2;  % Faktor eksponencijalnog prilagođavanja za brži prelaz
beta = 1.5; % Damping faktor za kontrolu intenziteta

%% Prelaz od R ka G (Transition R → G)
R_to_G = zeros(cols, 3); % Svaka tačka sadrži [R, G, B]
R_to_G(:,1) = exp(-alpha * (1 - lorenzCurve)); % Eksponencijalno smanjenje crvene komponente
R_to_G(:,2) = 1 - exp(-beta * lorenzCurve);    % Eksponencijalno povećanje zelene komponente
R_to_G(:,3) = 0; % Plava komponenta ostaje 0

%% Prelaz od G ka B (Transition G → B)
G_to_B = zeros(cols, 3); % Svaka tačka sadrži [R, G, B]
G_to_B(:,1) = 0; % Crvena komponenta ostaje 0
G_to_B(:,2) = exp(-alpha * (1 - lorenzCurve)); % Eksponencijalno smanjenje zelene komponente
G_to_B(:,3) = 1 - exp(-beta * lorenzCurve);    % Eksponencijalno povećanje plave komponente

%% Prelaz od B ka R (Transition B → R)
B_to_R = zeros(cols, 3); % Svaka tačka sadrži [R, G, B]
B_to_R(:,1) = 1 - exp(-beta * lorenzCurve);    % Eksponencijalno povećanje crvene komponente
B_to_R(:,2) = 0; % Zelena komponenta ostaje 0
B_to_R(:,3) = exp(-alpha * (1 - lorenzCurve)); % Eksponencijalno smanjenje plave komponente

%% Crtanje 3D RGB prelaza
figure;

% R → G prelaz
plot3(R_to_G(:,1), R_to_G(:,2), R_to_G(:,3), 'r-', 'LineWidth', 2); hold on;

% G → B prelaz
plot3(G_to_B(:,1), G_to_B(:,2), G_to_B(:,3), 'g-', 'LineWidth', 2);

% B → R prelaz
plot3(B_to_R(:,1), B_to_R(:,2), B_to_R(:,3), 'b-', 'LineWidth', 2);

% Podešavanje 3D grafika
grid on;
xlabel('Crvena (R) / Red (R)');
ylabel('Zelena (G) / Green (G)');
zlabel('Plava (B) / Blue (B)');
title('3D RGB prelazi sa Lorencovom krivom i diferencijalnim jednačinama');
legend({'R → G', 'G → B', 'B → R'}, 'Location', 'Best');
view(3); % 3D prikaz