clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script per risolvere numericamente il campo di temperatura 
% esercizio numero 1.53 
% metodi: Eulero, Laasonen, Crank-Nicholson, confronto
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1- ask for the parameters
prompt={'Numero di nodi',... %1
    'Diffusività termica [m^2/s]',... %2
    'Conducibilità termica [W/mK]',... %3
    'Lunghezza aletta [m]',... %4
        'Altezza saldatore [m]',... %5
        'Spessore saldatore [m]',... %6
        'Temperatura iniziale corpo [°C]',... %7
        'Temperatura ambiente [°C]',... %8
        'Coefficiente di scambio convettivo [W/m^2K]',... %9
        'Generazione interna [W]',... %10
        'Potenza scambiata dalla base, uscente positiva [W]',...%11
        'Metodo: [1] Eulero, [2] Laasonen, [3] Crank-Nicholson, [4] Confronto'}; %12
name='Input dati';
numlines=1;
defaultanswer={'100','1.53e-5','55','0.100','0.003','0.003','25','25','15','50','10','4'};
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';
answerf=inputdlg(prompt,name,numlines,defaultanswer,options);

% Estrazione dati dall'input
N = str2double(cell2mat(answerf(1)));
alpha = str2double(cell2mat(answerf(2)));
k = str2double(cell2mat(answerf(3)));
s = str2double(cell2mat(answerf(4)));
b = str2double(cell2mat(answerf(6)));
h = str2double(cell2mat(answerf(5)));
T_in = str2double(cell2mat(answerf(7)));
T_amb = str2double(cell2mat(answerf(8)));
h_conv = str2double(cell2mat(answerf(9)));
U_gen = str2double(cell2mat(answerf(10)));
Q_base = str2double(cell2mat(answerf(11)));
metodo = str2double(cell2mat(answerf(12)));

% CALCOLO DEL DX
dx = s / (N - 1);

%% SOLUZIONI

%% [1] Metodo di Eulero

if metodo == 1 % Eulero
    % Limiti di stabilità
    dt_max1 = dx * dx / (2 * alpha);
    dt_max2 = (dx * dx / (2 * alpha)) / (1 + h_conv * dx / k);
    dt_max=round(min(dt_max1,dt_max2), 3);
    dt_max_str = num2str(dt_max);
    % Richiesta dt
    prompt = {['Passo temporale? (Necessariamente minore di ', dt_max_str, ') [s]'], ... %1
        'Tempo di fine simulazione? [s]'}; %2
    name = 'Input dati';
    numlines = 1;
    defaultanswer = {'0.01', '110'}; 
    options.Resize = 'on';
    options.WindowStyle = 'normal';
    options.Interpreter = 'tex';
    answerf2 = inputdlg(prompt, name, numlines, defaultanswer, options);    
    dt = str2double(cell2mat(answerf2(1)));
    t_end = str2double(cell2mat(answerf2(2)));

    % Variabili di calcolo
    dFo = alpha * dt / (dx^2);
    U3 = U_gen/(b * h * s);
    rho_c = k / alpha;
    Bi = h_conv * dx / k;
    M = floor(t_end / dt); %suddivisioni temporali
    % Imposizione condizioni iniziali
    T = zeros(N, M+1);
    for I = 1:N
        T(I, 1) = T_in;
    end
    for j = 1:M
        TIME = j * dt;

        % PRIMO NODO
        T(1, j+1) = T(1, j) + 2 * dFo * (T(2, j) - T(1, j)) + (4 * h_conv * dt * ((T_amb - T(1, j))) / (rho_c * b)) + (U3 * dt / rho_c) - 2* dt * (Q_base * (1 - exp(-0.050 * dt * j)) / (rho_c * dx * b * h));

        % NODI INTERNI
        for I = 2:(N - 1)
            T(I, j+1) = dFo * (T(I - 1, j) + T(I + 1, j)) + (1 - 2 * dFo) * T(I, j) + (U3 * dt / rho_c) + ((4 * h_conv * dt * (T_amb - T(I, j))) / (rho_c * b));
        end

        % ULTIMO NODO
        T(N, j+1) = T(N, j) + 2 * dFo * (T(N - 1, j) - T(N, j)) + (4 * h_conv * dt * ((T_amb - T(N, j)))/(rho_c * b)) + ((2 * h_conv * dt * (T_amb - T(N, j))) / (rho_c * dx)) + (U3 * dt / (rho_c));
    end

    % Plot delle temperature
    num_nodi_plot = 5;
    nodi_plot = [1, floor(N/4), floor(N/2), floor(3*N/4), N];
    if floor(N/4) == 1
        num_nodi_plot = 4;
        nodi_plot = [floor(N/4), floor(N/2), floor(3*N/4), N];
    end

    K = T(nodi_plot, :);
    t = linspace(0, t_end, M+1);

    figure
    plot(t, K)
    hold on
    yline(1200, '-', '1200\circC')

    % Aggiunta della legenda
    labels = cell(num_nodi_plot, 1);
    for i = 1:num_nodi_plot
        labels{i} = ['Nodo ' num2str(nodi_plot(i))];
    end
    labels{end+1} = 'T desiderata [°C]';
    legend(labels, 'Location', 'best')

    % Grafico live

    % Creazione della figura per l'animazione
    x = linspace(0, s, N);
    figure;
    h_plot_eulero = plot(x, T(:, 1));
    xlabel('Posizione x lungo il saldatore [m]');
    ylabel('T [°C]');
    title('T(x), Eulero');

    % Imposta i limiti dell'asse Y in modo fisso
    y_min_eulero = min(T(:)); % Valore minimo della temperatura
    y_max_eulero = max(T(:)); % Valore massimo della temperatura
    ylim([y_min_eulero, 1.1*y_max_eulero]);

    for j = 2:M+1
        % Aggiornamento del grafico animato
        set(h_plot_eulero, 'YData', T(:, j));
        title(['Tempo = ' num2str((j-1)*dt) ' s']);

        % Aggiornamento della visualizzazione
        drawnow;

        % Pausa per controllare la velocità dell'animazione
        pause(0.005); % Modifica il valore della pausa se desideri cambiare la velocità dell'animazione
    end

    % Mostra il grafico completo alla fine dell'animazione
    figure;
    plot(x, T(:, end));
    xlabel('Posizione lungo il saldatore [m]');
    ylabel('T [°C]');
    title('T(x) a \theta_{end}, Eulero');

%% Metodo di Lasoonen [2]
elseif metodo==2 %Laasonen
    prompt={'Passo temporale? [s]',... %1
        'Tempo di fine simulazione? [s]'}; %2
    name='Input dati';
    numlines=1;
    defaultanswer={'0.01','110'};
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
    answerf2=inputdlg(prompt,name,numlines,defaultanswer,options);
    dt = str2double(cell2mat(answerf2(1)));
    t_end = str2double(cell2mat(answerf2(2)));

    % Variabili di calcolo
    dFo = alpha * dt / (dx^2);
    U3 = U_gen/(b * h * s);
    rho_c = k / alpha;
    Bi = h_conv * dx / k;
    M = floor(t_end / dt); %suddivisioni temporali
    % Coefficienti metodo
    l_1 = -alpha * dt / dx^2;
    l_3 = -4 * l_1 * Bi * dx / b;
    l_2 = 1 - 2 * l_1 + l_3;
    l_4 = alpha * dt / k;
    l_5 = 2 * alpha * Bi * dt; 

    % Condizioni iniziali
    T_0 = T_in * ones(N, 1);
    % Risoluzione del sistema lineare per ogni passo temporale
    T = zeros(N, M+1);
    T(:, 1) = T_0;
    z = zeros(N, M+1); % Inizializzazione della variabile z
    for j = 2:M+1
    A = diag(l_2 * ones(N, 1));
        for i = 1:N-1
            A(i, i+1) = l_1;
        end
        for i = 1:N-1
            A(i+1, i) = l_1;
        end
    A(1, 2) = 2 * l_1;
    A(N, N-1) = 2 * l_1;
    A(N, N) = l_2 + l_5;
    z(1:N, j) = T(1:N, j-1) + l_3 * T_amb + l_4 * U3;
    z(1, j) = T(1, j-1) + l_3 * T_amb + l_4 * U3 - ((2 * dx * dFo) / (k * b * h)) * (Q_base * (1 - exp(-0.05 * dt * (j-1))));
    z(N, j) = T(N, j-1) + (l_3 + l_5) * T_amb + l_4 * U3;
    T(1:N, j) = thomas(A, z(:, j));
    end

    % PLOT TEMPORALE
    num_nodi_plot = 5;
    nodi_plot = [1, floor(N/4), floor(N/2), floor(3*N/4), N];
    if floor(N/4) == 1
        num_nodi_plot = 4;
        nodi_plot = [floor(N/4), floor(N/2), floor(3*N/4), N];
    end

    
    K = T(nodi_plot, :);
    t = linspace(0, t_end, M+1);
    
    figure
    plot(t, K)
    hold on
    yline(1200, '-', '1200\circC')
    
    % Aggiunta della legenda
    labels = cell(num_nodi_plot, 1);
    for i = 1:num_nodi_plot
        labels{i} = ['Nodo ' num2str(nodi_plot(i))];
    end
    labels{end+1} = 'T desiderata [°C]';
    legend(labels, 'Location', 'best')
    
    % Grafico live
    % Creazione della figura per l'animazione
    x = linspace(0, s, N);
    figure;
    h_plot = plot(x, T(:, 1));
    xlabel('Posizione x lungo il saldatore [m]');
    ylabel('T [°C]');
    title('T(x), Laasonen');
    
    % Imposta i limiti dell'asse Y in modo fisso
    y_min = min(T(:)); % Valore minimo della temperatura
    y_max = max(T(:)); % Valore massimo della temperatura
    ylim([y_min, 1.1*y_max]);
    
    for j = 2:M+1
        % Aggiornamento del grafico animato
        set(h_plot, 'YData', T(:, j));
        title(['Tempo = ' num2str((j-1)*dt) ' s']);
        
        % Aggiornamento della visualizzazione
        drawnow;
        
        % Pausa per controllare la velocità dell'animazione
        pause(0.005); 
    end
    
    % Mostra il grafico completo alla fine dell'animazione
    figure;
    plot(x, T(:, end));
    xlabel('Posizione x lungo il saldatore [m]');
    ylabel('T [°C]');
    title('T(x) a \theta_{end}, Laasonen');

%% Crank-Nicholson [3] 
elseif metodo==3 %Crank-Nicholson    
    U3 = U_gen/(b * h * s);
    rho_c = k / alpha;
    Bi = h_conv * dx / k;
    % Limiti di stabilità
    dt_max1 = ((1)/(((alpha)/(dx^2))+((2*alpha)/(b*dx))));
    dt_max2 = ((1)/(((2*h_conv)/(rho_c*b))+((alpha)/(dx^2))));
    dt_max = round(min(dt_max1, dt_max2), 3);
    dt_max_str = num2str(dt_max);
    prompt = {['Passo temporale? (Necessariamente minore di ', dt_max_str, ') [s]'], ... %1
        'Tempo di fine simulazione? [s]'}; %2
    
    name = 'Input dati';
    numlines = 1;
    defaultanswer = {'0.01', '110'};
    options.Resize = 'on';
    options.WindowStyle = 'normal';
    options.Interpreter = 'tex';
    answerf2 = inputdlg(prompt, name, numlines, defaultanswer, options);
    dt = str2double(cell2mat(answerf2(1)));
    t_end = str2double(cell2mat(answerf2(2)));
   
    % Variabili di calcolo
    dFo = alpha * dt / (dx^2);
    M = floor(t_end / dt); %suddivisioni temporali

    % Coefficienti
    c_1 = -dFo;
    c_2 = 2*dFo*Bi*dx/b;
    c_3 = (1-c_1+c_2);
    c_4 = alpha*Bi*dt;

    % Creazione della matrice di coefficienti per il metodo di Thomas
    A = diag(c_3 * ones(N, 1));
    for i = 1:N-1
        A(i, i+1) = c_1/2;
    end
    for i = 1:N-1
        A(i+1, i) = c_1/2;
    end
    A(1, 2) = c_1;
    A(N, N-1) = c_1;
    A(N, N) = (c_3+c_4);

    % Condizioni iniziali
    T_0 = T_in * ones(N, 1);
    % Risoluzione del sistema lineare per ogni passo temporale
    T = zeros(N, M+1);
    T(:, 1) = T_0;
    z = zeros(N, M+1); % Inizializzazione della variabile z
    for j = 2:M+1
        z(2:N-1, j) = -c_1/2*(T(1:N-2, j-1)+T(3:N, j-1)) + (1+c_1-c_2)*T(2:N-1, j-1) - 1*c_2*T_amb + (alpha*dt/k)*U3;
        z(1, j) = (1+c_1-c_2)*T(1, j-1) - c_1*T(2, j-1) + 2*c_2*T_amb + (alpha*dt/k)*U3 - ((2*dx*dFo)/(b*h*k)) * 0.5 *((Q_base * (1 - exp(-0.05 * dt * (j)))) + (Q_base * (1 - exp(-0.05 * dt * (j-1)))));
        z(N, j) = -c_1*T(N-1, j-1) + (1+c_1-c_2-c_4)*T(N, j-1) + (2*(c_2+c_4))*T_amb + (alpha*dt/k)*U3;
        T(1:N, j) = thomas(A, z(:, j));
    end

    % PLOT TEMPORALE
    num_nodi_plot = 5;
    nodi_plot = [1, floor(N/4), floor(N/2), floor(3*N/4), N];
    if floor(N/4) == 1
        num_nodi_plot = 4;
        nodi_plot = [floor(N/4), floor(N/2), floor(3*N/4), N];
    end
   
    K = T(nodi_plot, :);
    t = linspace(0, t_end, M+1);
    
    figure
    plot(t, K)
    hold on
    yline(1200, '-', '1200\circC')
    
    % Aggiunta della legenda
    labels = cell(num_nodi_plot, 1);
    for i = 1:num_nodi_plot
        labels{i} = ['Nodo ' num2str(nodi_plot(i))];
    end
    labels{end+1} = 'T desiderata [°C]';
    legend(labels, 'Location', 'best')
    
    % Trova gli indici dei punti di intersezione
    intersezioni = find(K == 1200);
    
    % Aggiungi i punti di intersezione al grafico
    plot(t(intersezioni), K(intersezioni), 'ro', 'MarkerSize', 8, 'LineWidth', 2)
    
    hold off
    xlabel('\theta [s]')
    ylabel('T [°C]')
    title('T(\theta), Crank-Nicholson')
    
    % Grafico live
    
    % Creazione della figura per l'animazione
    x = linspace(0, s, N);
    figure;
    h_plot = plot(x, T(:, 1));
    xlabel('Posizione x lungo il saldatore [m]');
    ylabel('T [°C]');
    title('T(x), Crank-Nicholson');
    
    % Imposta i limiti dell'asse Y in modo fisso
    y_min = min(T(:)); % Valore minimo della temperatura
    y_max = max(T(:)); % Valore massimo della temperatura
    ylim([y_min, y_max + (0.1*y_max)]);
    
    for j = 2:M+1
        % Aggiornamento del grafico animato
        set(h_plot, 'YData', T(:, j));
        title(['Tempo = ' num2str((j-1)*dt) ' s']);
        
        % Aggiornamento della visualizzazione
        drawnow;
        
        % Pausa per controllare la velocità dell'animazione
        pause(0.005); % Modifica il valore della pausa se desideri cambiare la velocità dell'animazione
    end
    
    % Mostra il grafico completo alla fine dell'animazione
    figure;
    plot(x, T(:, end));
    xlabel('Posizione x lungo il saldatore [m]');
    ylabel('T [°C]');
    title('T(x) a \theta_{end}, Crank-Nicholson');
elseif metodo==4
    U3 = U_gen/(b * h * s);
    rho_c = k / alpha;
    Bi = h_conv * dx / k;
    %stabilità
    dt_max1_cn = ((1)/(((alpha)/(dx^2))+((2*alpha)/(b*dx))));
    dt_max2_cn = ((1)/(((2*h_conv)/(rho_c*b))+((alpha)/(dx^2))));
    dt_max1_eu = dx * dx / (2 * alpha);
    dt_max2_eu = (dx * dx / (2 * alpha)) / (1 + h_conv * dx / k);
    dt_max_cn = min(dt_max1_cn,dt_max2_cn);
    dt_max_eu = min(dt_max1_eu,dt_max2_eu);
    dt_max = round(min(dt_max_cn, dt_max_eu), 3);
    dt_max_str = num2str(dt_max);
    % Richiesta dt
    prompt = {['Passo temporale? (Necessariamente minore di ', dt_max_str, ') [s]'], ... %1
        'Tempo di fine simulazione? [s]'}; %2
    name = 'Input dati';
    numlines = 1;
    defaultanswer = {'0.01', '110'}; % Assicurati che ci siano al massimo due elementi
    options.Resize = 'on';
    options.WindowStyle = 'normal';
    options.Interpreter = 'tex';
    answerf2 = inputdlg(prompt, name, numlines, defaultanswer, options);    
    dt = str2double(cell2mat(answerf2(1)));
    t_end = str2double(cell2mat(answerf2(2)));

    %1
    % Variabili di calcolo
    dFo = alpha * dt / (dx^2);
    M = floor(t_end / dt); %suddivisioni temporali
    % Imposizione condizioni iniziali
    T_eu = zeros(N, M+1);
    for I = 1:N
        T_eu(I, 1) = T_in;
    end
    for j = 1:M
        TIME = j * dt;

        % PRIMO NODO
        T_eu(1, j+1) = T_eu(1, j) + 2 * dFo * (T_eu(2, j) - T_eu(1, j)) + (4 * h_conv * dt * ((T_amb - T_eu(1, j))) / (rho_c * b)) + (U3 * dt / rho_c) - 2*dt * (Q_base * (1 - exp(-0.050 * dt * j)) / (rho_c * dx * b * h));

        % NODI INTERNI
        for I = 2:(N - 1)
            T_eu(I, j+1) = dFo * (T_eu(I - 1, j) + T_eu(I + 1, j)) + (1 - 2 * dFo) * T_eu(I, j) + (U3 * dt / rho_c) + ((4 * h_conv * dt * (T_amb - T_eu(I, j))) / (rho_c * b));
        end

        % ULTIMO NODO
        T_eu(N, j+1) = T_eu(N, j) + 2 * dFo * (T_eu(N - 1, j) - T_eu(N, j)) + (4 * h_conv * dt * ((T_amb - T_eu(N, j)))/(rho_c * b)) + ((2 * h_conv * dt * (T_amb - T_eu(N, j))) / (rho_c * dx)) + (U3 * dt / (rho_c));
    end
    %2
    % Coefficienti metodo
    l_1 = -alpha * dt / dx^2;
    l_3 = -4 * l_1 * Bi * dx / b;
    l_2 = 1 - 2 * l_1 + l_3;
    l_4 = alpha * dt / k;
    l_5 = 2 * alpha * Bi * dt; 
    % Creazione della matrice di coefficienti per il metodo di Thomas
    A = diag(l_2 * ones(N, 1));
    for i = 1:N-1
        A(i, i+1) = l_1;
    end
    for i = 1:N-1
        A(i+1, i) = l_1;
    end
    A(1, 2) = 2 * l_1;
    A(N, N-1) = 2 * l_1;
    A(N, N) = l_2 + l_5;
    % Condizioni iniziali
    T_las0 = T_in * ones(N, 1);
    % Risoluzione del sistema lineare per ogni passo temporale
    T_las = zeros(N, M+1);
    T_las(:, 1) = T_las0;
    z = zeros(N, M+1); % Inizializzazione della variabile z
    for j = 2:M+1
        z(1:N, j) = T_las(1:N, j-1) + l_3 * T_amb + l_4 * U3;
        z(1, j) = T_las(1, j-1) + l_3 * T_amb + l_4 * U3 - ((2 * dx * dFo) / (k * b * h)) * (Q_base * (1 - exp(-0.05 * dt * (j-1))));
        z(N, j) = T_las(N, j-1) + (l_3 + l_5) * T_amb + l_4 * U3;
        T_las(1:N, j) = thomas(A, z(:, j));
    end
    %3
    % Coefficienti
    c_1 = -dFo;
    c_2 = 2*dFo*Bi*dx/b;
    c_3 = (1-c_1+c_2);
    c_4 = alpha*Bi*dt;

    % Creazione della matrice di coefficienti per il metodo di Thomas
    A = diag(c_3 * ones(N, 1));
    for i = 1:N-1
        A(i, i+1) = c_1/2;
    end
    for i = 1:N-1
        A(i+1, i) = c_1/2;
    end
    A(1, 2) = c_1;
    A(N, N-1) = c_1;
    A(N, N) = (c_3+c_4);

    % Condizioni iniziali
    T_cn0 = T_in * ones(N, 1);
    % Risoluzione del sistema lineare per ogni passo temporale
    T_cn = zeros(N, M+1);
    T_cn(:, 1) = T_cn0;
    z = zeros(N, M+1); % Inizializzazione della variabile z
    for j = 2:M+1
        z(2:N-1, j) = -c_1/2*(T_cn(1:N-2, j-1)+T_cn(3:N, j-1)) + (1+c_1-c_2)*T_cn(2:N-1, j-1) - 1*c_2*T_amb + (alpha*dt/k)*U3;
        z(1, j) = (1+c_1-c_2)*T_cn(1, j-1) - c_1*T_cn(2, j-1) + 2*c_2*T_amb + (alpha*dt/k)*U3 - ((2*dx*dFo)/(b*h*k)) * 0.5 *((Q_base * (1 - exp(-0.05 * dt * (j)))) + (Q_base * (1 - exp(-0.05 * dt * (j-1)))));
        z(N, j) = -c_1*T_cn(N-1, j-1) + (1+c_1-c_2-c_4)*T_cn(N, j-1) + (2*(c_2+c_4))*T_amb + (alpha*dt/k)*U3;
        T_cn(1:N, j) = thomas(A, z(:, j));
    end    
    % Plot delle temperature
    figure;
           
    x = linspace(0, t_end, M+1);
    plot(x, T_eu(end, :));
    hold on
    plot(x, T_las(end, :));
    hold on
    plot(x, T_cn(end, :));
    hold on
    yline(1200, '-', '1200\degC')
    hold off

    % Imposta i limiti dell'asse Y in modo fisso
    y_min = min(T_las(:)); % Valore minimo della temperatura
    y_max = max(T_las(:)); % Valore massimo della temperatura
    ylim([y_min, y_max + (0.1*y_max)]);

    legend({'Eulero', 'Laasonen', 'Crank-Nicolson', 'T desiderata [°C]'}, 'Location', 'southwest')
    xlabel('Tempo \theta [s]');
    ylabel('T [°C]');
    title('T(\theta) per i tre metodi');

else
    disp('Inserire un numero valido.')
end









