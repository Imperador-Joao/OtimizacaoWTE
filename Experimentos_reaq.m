clc
clear all
close all

addpath("Funções")

tic

%Tabelas

LiqeVap = readtable("../Tabelas/Liq e vapor saturado.xlsx", "TextType", "string");
Liqcomp = readtable("../Tabelas/Liq subresfriado.xlsx", "TextType", "string");
VapSuperAq = readtable("../Tabelas/Vapor superaquecido.xlsx", "TextType", "string");
Sat = readtable("../Tabelas/Curva de saturação.xlsx", "TextType", "string");

%Dados do ciclo

T_comp = Liqcomp.Temperature_C_;
P_comp = Liqcomp.Pressure_bar_;
h_comp = Liqcomp.Enthalpy_kJ_kg_;
s_comp = Liqcomp.Entropy_J_g_K_;

T_sat = LiqeVap.Temperature_C_;
P_sat = LiqeVap.Pressure_bar_;
h_sat = LiqeVap.Enthalpy_kJ_kg_;
s_sat = LiqeVap.Entropy_J_g_K_;
estado_fisico = LiqeVap.Phase;

T_superaq = VapSuperAq.Temperature_C_;
P_superaq = VapSuperAq.Pressure_bar_;
h_superaq = VapSuperAq.Enthalpy_kJ_kg_;
s_superaq = VapSuperAq.Entropy_J_g_K_;

T_sat = T_sat(estado_fisico == "liquid");
P_sat = P_sat(estado_fisico == "liquid");

h_sat_liq = h_sat(estado_fisico == "liquid");
h_sat_vap = h_sat(estado_fisico == "vapor");

s_sat_liq = s_sat(estado_fisico == "liquid");
s_sat_vap = s_sat(estado_fisico == "vapor");

Pmin = 0.12;
Pmax = 100;

Ef_isoentropica = 0.9;

%Otimizar

Matriz_Entropias_Temperaturas = [];
Matriz_Pressoes_Temp_reaq_Eficiencia_x14 = [];
i = 1;

%P_inter1 = 5;
%P_inter2 = 10;
%P_reaq = 15;
%P_inter3 = 20;
%T_reaq = 490;

for P_inter1 = 5:5:95                                                      % Pressão da extração 1
    for P_inter2 = 10:5:95
        
        if P_inter2 <= P_inter1                                            %P_inter2 deve ser maior que P_inter1
            continue
        end
        
        for P_reaq = 15:5:95
            
            if P_reaq <= P_inter2
                continue
            end
            
            for P_inter3 = 20:5:95
                
                if P_inter3 <= P_reaq
                    continue
                end
                
                % P_inter1 < P_inter2 < P_reaq < P_inter3
                
                for T_reaq = 260:20:500

                    T_comp_P_inter1 = T_comp(P_comp == P_inter1);
                    h_comp_P_inter1 = h_comp(P_comp == P_inter1);
                    s_comp_P_inter1 = s_comp(P_comp == P_inter1);

                    T_comp_Pmax = T_comp(P_comp == Pmax);
                    h_comp_Pmax = h_comp(P_comp == Pmax);
                    s_comp_Pmax = s_comp(P_comp == Pmax);

                    T_superaq_Pmin = T_superaq(P_superaq == Pmin);
                    h_superaq_Pmin = h_superaq(P_superaq == Pmin);
                    s_superaq_Pmin = s_superaq(P_superaq == Pmin);


                    T_superaq_P_inter1 = T_superaq(P_superaq == P_inter1);
                    h_superaq_P_inter1 = h_superaq(P_superaq == P_inter1);
                    s_superaq_P_inter1 = s_superaq(P_superaq == P_inter1);

                    T_superaq_P_inter2 = T_superaq(P_superaq == P_inter2);
                    h_superaq_P_inter2 = h_superaq(P_superaq == P_inter2);
                    s_superaq_P_inter2 = s_superaq(P_superaq == P_inter2);

                    T_superaq_P_reaq = T_superaq(P_superaq == P_reaq);
                    h_superaq_P_reaq = h_superaq(P_superaq == P_reaq);
                    s_superaq_P_reaq = s_superaq(P_superaq == P_reaq);

                    T_superaq_P_inter3 = T_superaq(P_superaq == P_inter3);
                    h_superaq_P_inter3 = h_superaq(P_superaq == P_inter3);
                    s_superaq_P_inter3 = s_superaq(P_superaq == P_inter3);


                    T_superaq_Pmax = T_superaq(P_superaq == Pmax);
                    h_superaq_Pmax = h_superaq(P_superaq == Pmax);
                    s_superaq_Pmax = s_superaq(P_superaq == Pmax);

                    s_sat_Pmin = s_superaq_Pmin(1);
                    s_sat_P_inter1 = s_superaq_P_inter1(1);
                    s_sat_P_inter2 = s_superaq_P_inter2(1);
                    s_sat_P_reaq = s_superaq_P_reaq(1);
                    s_sat_P_inter3 = s_superaq_P_inter3(1);

                    h_sat_Pmin = h_superaq_Pmin(1);
                    h_sat_P_inter1 = h_superaq_P_inter1(1);
                    h_sat_P_inter2 = h_superaq_P_inter2(1);
                    h_sat_P_reaq = h_superaq_P_reaq(1);
                    h_sat_P_inter3 = h_superaq_P_inter3(1);

                    %Ponto 1

                    P1 = Pmin;
                    T1 = T_sat(1);
                    h1 = h_sat_liq(1);
                    s1 = s_sat_liq(1);

                    %Ponto 2s

                    P2s = P_inter1;
                    s2s = s1;
                    T2s = Interpolar(s_comp_P_inter1,T_comp_P_inter1,s2s);
                    h2s = Interpolar(s_comp_P_inter1,h_comp_P_inter1,s2s);

                    %Ponto 2

                    P2 = P2s;
                    h2 = h1+(h2s-h1)/Ef_isoentropica;
                    T2 = Interpolar(h_comp_P_inter1,T_comp_P_inter1,h2);
                    s2 = Interpolar(h_comp_P_inter1,s_comp_P_inter1,h2);


                    %Ponto 3

                    P3 = P2;
                    T3 = Interpolar(P_sat,T_sat,P3);
                    h3 = Interpolar(P_sat,h_sat_liq,P3);
                    s3 = Interpolar(P_sat,s_sat_liq,P3);

                    %Ponto 4s

                    P4s = Pmax;
                    s4s = s3;
                    T4s = Interpolar(s_comp_Pmax,T_comp_Pmax,s4s);
                    h4s = Interpolar(s_comp_Pmax,h_comp_Pmax,s4s);

                    %Ponto 4

                    P4 = P4s;
                    h4 = h3+(h4s-h3)/Ef_isoentropica;
                    T4 = Interpolar(h_comp_Pmax,T_comp_Pmax,h4);
                    s4 = Interpolar(h_comp_Pmax,s_comp_Pmax,h4);

                    %Ponto 5

                    P5 = P4;
                    T5 = Interpolar(P_sat,T_sat,P_inter2);
                    h5 = Interpolar(T_comp_Pmax,h_comp_Pmax,T5);
                    s5 = Interpolar(T_comp_Pmax,s_comp_Pmax,T5);

                    % Ponto 6

                    P6 = P5;
                    T6 = Interpolar(P_sat,T_sat,P_inter3);
                    h6 = Interpolar(T_comp_Pmax,h_comp_Pmax,T6);
                    s6 = Interpolar(T_comp_Pmax,s_comp_Pmax,T6);

                    %Ponto 7

                    P7 = P4;
                    T7 = T_sat(end);
                    h7 = h_sat_vap(end);
                    s7 = s_sat_vap(end);

                    %Ponto 8

                    P8 = P7;
                    T8 = 500;
                    h8 = Interpolar(T_superaq_Pmax,h_superaq_Pmax,T8);
                    s8 = Interpolar(T_superaq_Pmax,s_superaq_Pmax,T8);

                    %Ponto 9s

                    P9s = P_inter3;
                    s9s = s8;

                    if s9s > s_sat_P_inter3                                                 % Verificar se o ponto está fora da curva de saturação
                        h9s = Interpolar(s_superaq_P_inter3,h_superaq_P_inter3,s9s);
                        T9s = Interpolar(s_superaq_P_inter3,T_superaq_P_inter3,s9s);
                    else                                                                    %O ponto estando dentro da curva de saturação, o seu título será calculado
                        T9s = Interpolar(P_sat,T_sat,P9s);
                        s9s_liq = Interpolar(P_sat,s_sat_liq,P9s);
                        s9s_vap = Interpolar(P_sat,s_sat_vap,P9s);
                        x9s = (s9s-s9s_liq)/(s9s_vap-s9s_liq);

                        h9s_liq = Interpolar(P_sat,h_sat_liq,P9s);
                        h9s_vap = Interpolar(P_sat,h_sat_vap,P9s);
                        h9s = (1-x9s)*h9s_liq + x9s*h9s_vap;
                    end

                   %Ponto 9

                   P9 = P9s;
                   h9 = h8 - (h8 - h9s)*Ef_isoentropica;

                   if h9 > h_sat_P_inter3
                       s9 = Interpolar(h_superaq_P_inter3,s_superaq_P_inter3,h9);
                       T9 = Interpolar(h_superaq_P_inter3,T_superaq_P_inter3,h9);
                   else
                        T9 = Interpolar(P_sat,T_sat,P9);
                        h9_liq = Interpolar(P_sat,h_sat_liq,P9);
                        h9_vap = Interpolar(P_sat,h_sat_vap,P9);
                        x9 = (h9-h9_liq)/(h9_vap-h9_liq);

                        s9_liq = Interpolar(P_sat,s_sat_liq,P9);
                        s9_vap = Interpolar(P_sat,s_sat_vap,P9);
                        s9 = (1-x9)*s9_liq + x9*s9_vap;

                   end

                   %Ponto 10s

                   P10s = P_reaq;
                   s10s = s9;

                    if s10s > s_sat_P_reaq
                        h10s = Interpolar(s_superaq_P_reaq,h_superaq_P_reaq,s10s);
                        T10s = Interpolar(s_superaq_P_reaq,T_superaq_P_reaq,s10s);
                    else
                        T10s = Interpolar(P_sat,T_sat,P10s);
                        s10s_liq = Interpolar(P_sat,s_sat_liq,P10s);
                        s10s_vap = Interpolar(P_sat,s_sat_vap,P10s);
                        x10s = (s10s-s10s_liq)/(s10s_vap-s10s_liq);

                        h10s_liq = Interpolar(P_sat,h_sat_liq,P10s);
                        h10s_vap = Interpolar(P_sat,h_sat_vap,P10s);
                        h10s = (1-x10s)*h10s_liq + x10s*h10s_vap;
                    end

                   %Ponto 10

                   P10 = P10s;
                   h10 = h9 - (h9 - h10s)*Ef_isoentropica;

                   if h10 > h_sat_P_reaq
                       s10 = Interpolar(h_superaq_P_reaq,s_superaq_P_reaq,h10);
                       T10 = Interpolar(h_superaq_P_reaq,T_superaq_P_reaq,h10);
                   else
                        T10 = Interpolar(P_sat,T_sat,P10);
                        h10_liq = Interpolar(P_sat,h_sat_liq,P10);
                        h10_vap = Interpolar(P_sat,h_sat_vap,P10);
                        x10 = (h10-h10_liq)/(h10_vap-h10_liq);

                        s10_liq = Interpolar(P_sat,s_sat_liq,P10);
                        s10_vap = Interpolar(P_sat,s_sat_vap,P10);
                        s10 = (1-x10)*s10_liq + x10*s10_vap;

                   end
                   
                   if T10 >= T_reaq
                       continue
                   end

                   %Ponto 11

                   P11 = P10;
                   T11 = T_reaq;
                   h11 = Interpolar(T_superaq_P_reaq,h_superaq_P_reaq,T11);
                   s11 = Interpolar(T_superaq_P_reaq,s_superaq_P_reaq,T11);

                   %Ponto 12s

                    P12s = P_inter2;
                    s12s = s11;

                    if s12s > s_sat_P_inter2
                        h12s = Interpolar(s_superaq_P_inter2,h_superaq_P_inter2,s12s);
                        T12s = Interpolar(s_superaq_P_inter2,T_superaq_P_inter2,s12s);
                    else
                        T12s = Interpolar(P_sat,T_sat,P12s);
                        s12s_liq = Interpolar(P_sat,s_sat_liq,P12s);
                        s12s_vap = Interpolar(P_sat,s_sat_vap,P12s);
                        x12s = (s12s-s12s_liq)/(s12s_vap-s12s_liq);

                        h12s_liq = Interpolar(P_sat,h_sat_liq,P12s);
                        h12s_vap = Interpolar(P_sat,h_sat_vap,P12s);
                        h12s = (1-x12s)*h12s_liq + x12s*h12s_vap;
                    end

                   %Ponto 12

                   P12 = P12s;
                   h12 = h11 - (h11 - h12s)*Ef_isoentropica;

                   if h12 > h_sat_P_inter2
                       s12 = Interpolar(h_superaq_P_inter2,s_superaq_P_inter2,h12);
                       T12 = Interpolar(h_superaq_P_inter2,T_superaq_P_inter2,h12);
                   else
                        T12 = Interpolar(P_sat,T_sat,P12);
                        h12_liq = Interpolar(P_sat,h_sat_liq,P12);
                        h12_vap = Interpolar(P_sat,h_sat_vap,P12);
                        x12 = (h12-h12_liq)/(h12_vap-h12_liq);

                        s12_liq = Interpolar(P_sat,s_sat_liq,P12);
                        s12_vap = Interpolar(P_sat,s_sat_vap,P12);
                        s12 = (1-x12)*s12_liq + x12*s12_vap;

                   end

                   %Ponto 13s

                    P13s = P_inter1;
                    s13s = s12;

                    if s13s > s_sat_P_inter1
                        h13s = Interpolar(s_superaq_P_inter1,h_superaq_P_inter1,s13s);
                        T13s = Interpolar(s_superaq_P_inter1,T_superaq_P_inter1,s13s);
                    else
                        T13s = Interpolar(P_sat,T_sat,P13s);
                        s13s_liq = Interpolar(P_sat,s_sat_liq,P13s);
                        s13s_vap = Interpolar(P_sat,s_sat_vap,P13s);
                        x13s = (s13s-s13s_liq)/(s13s_vap-s13s_liq);

                        h13s_liq = Interpolar(P_sat,h_sat_liq,P13s);
                        h13s_vap = Interpolar(P_sat,h_sat_vap,P13s);
                        h13s = (1-x13s)*h13s_liq + x13s*h13s_vap;
                    end

                   %Ponto 13

                   P13 = P13s;
                   h13 = h12 - (h12 - h13s)*Ef_isoentropica;

                   if h13 > h_sat_P_inter1
                       s13 = Interpolar(h_superaq_P_inter1,s_superaq_P_inter1,h13);
                       T13 = Interpolar(h_superaq_P_inter1,T_superaq_P_inter1,h13);
                   else
                        T13 = Interpolar(P_sat,T_sat,P13);
                        h13_liq = Interpolar(P_sat,h_sat_liq,P13);
                        h13_vap = Interpolar(P_sat,h_sat_vap,P13);
                        x13 = (h13-h13_liq)/(h13_vap-h13_liq);

                        s13_liq = Interpolar(P_sat,s_sat_liq,P13);
                        s13_vap = Interpolar(P_sat,s_sat_vap,P13);
                        s13 = (1-x13)*s13_liq + x13*s13_vap;

                   end           

                   %Ponto 14s

                    P14s = Pmin;
                    s14s = s13;

                    if s14s > s_sat_Pmin
                        h14s = Interpolar(s_superaq_Pmin,h_superaq_Pmin,s14s);
                        T14s = Interpolar(s_superaq_Pmin,T_superaq_Pmin,s14s);
                    else
                        T14s = T_sat(1);
                        s14s_liq = s_sat_liq(1);
                        s14s_vap = s_sat_vap(1);
                        x14s = (s14s-s14s_liq)/(s14s_vap-s14s_liq);

                        h14s_liq = h_sat_liq(1);
                        h14s_vap = h_sat_vap(1);
                        h14s = (1-x14s)*h14s_liq + x14s*h14s_vap;

                    end

                   %Ponto 14

                    P14 = P14s;
                    h14 = h13 - (h13 - h14s)*Ef_isoentropica;

                    if h14 > h_sat_Pmin
                       s14 = Interpolar(h_superaq_Pmin,s_superaq_Pmin,h14);
                       T14 = Interpolar(h_superaq_Pmin,T_superaq_Pmin,h14);
                    else
                        T14 = T_sat(1);
                        h14_liq = h_sat_liq(1);
                        h14_vap = h_sat_vap(1);
                        x14 = (h14-h14_liq)/(h14_vap-h14_liq);

                        s14_liq = s_sat_liq(1);
                        s14_vap = s_sat_vap(1);
                        s14 = (1-x14)*s14_liq + x14*s14_vap;

                    end

                   %Ponto 15

                    P15 = P9;
                    T15 = T6;
                    h15 = Interpolar(P_sat,h_sat_liq,P15);
                    s15 = Interpolar(P_sat,s_sat_liq,P15);

                   %Ponto 16

                    P16 = P12;
                    T16 = Interpolar(P_sat,T_sat,P16);
                    h16 = h15;

                    h16_liq = Interpolar(P_sat,h_sat_liq,P16);
                    h16_vap = Interpolar(P_sat,h_sat_vap,P16);
                    x16 = (h16-h16_liq)/(h16_vap-h16_liq);

                    s16_liq = Interpolar(P_sat,s_sat_liq,P16);
                    s16_vap = Interpolar(P_sat,s_sat_vap,P16);
                    s16 = (1-x16)*s16_liq+x16*s16_vap;

                   %Ponto 17

                    P17 = P16;
                    T17 = T5;
                    h17 = Interpolar(P_sat,h_sat_liq,P17);
                    s17 = Interpolar(P_sat,s_sat_liq,P17);

                   %Ponto 18

                    P18 = P13;
                    T18 = Interpolar(P_sat,T_sat,P18);
                    h18 = h17;

                    h18_liq = Interpolar(P_sat,h_sat_liq,P18);
                    h18_vap = Interpolar(P_sat,h_sat_vap,P18);
                    x18 = (h18-h18_liq)/(h18_vap-h18_liq);

                    s18_liq = Interpolar(P_sat,s_sat_liq,P18);
                    s18_vap = Interpolar(P_sat,s_sat_vap,P18);
                    s18 = (1-x18)*s18_liq+x18*s18_vap;

                    %Vazões mássicas

                    m = 1;
                    ma = m*(h6-h5)/(h9-h15);
                    mb = (m*(h5-h4)+ma*(h17-h16))/(h12-h17);
                    mc = (m*(h3-h2)+(ma+mb)*(h2-h18))/(h13-h2);
                    m_ma = m - ma;
                    m_ma_mb_mc = m - ma - mb - mc;
                    
                    if ma <= 0|| mb <= 0 || mc <= 0                        % As razões mássicas não podem ser negativas
                        continue
                    end

                    %Potências e eficiência

                    w12 = m_ma_mb_mc*(h2-h1);
                    w34 = m*(h4-h3);
                    w89 = m*(h8-h9);
                    w910 = m_ma*(h9-h10);
                    w1112 = m_ma*(h11-h12);
                    w1213 = (m_ma-mb)*(h12-h13);
                    w1314 = m_ma_mb_mc*(h13-h14);
                    Trabalho_turbina = w89+w910+w1112+w1213+w1314;
                    Trabalho_bombas = w12+w34;
                    Trabalho_liq = Trabalho_turbina-Trabalho_bombas;

                    q67 = m*(h7-h6);
                    q78 = m*(h8-h7);
                    q1011 = m_ma*(h11-h10);
                    Calor_fornecido = q67+q78+q1011;

                    Eficiencia = 100*Trabalho_liq/Calor_fornecido;

                    %Diagrama T-s

                    Entropias = [s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18];
                    Temperaturas = [T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17,T18];
                    
                    Matriz_Entropias_Temperaturas(:,:,i) = [Entropias;Temperaturas];
                    Matriz_Pressoes_Temp_reaq_Eficiencia_x14(i,:) = [P_inter1,P_inter2,P_reaq,P_inter3,T_reaq,Eficiencia,x14];
                    i = i+1;
                    
                end
            end
        end
    end
end

%Resultados

Eficiencias = Matriz_Pressoes_Temp_reaq_Eficiencia_x14(:,end-1);
[Maior_Ef,Pos] = max(Eficiencias)
Melhor_config = Matriz_Pressoes_Temp_reaq_Eficiencia_x14(Pos,:)
Melhores_entropias = Matriz_Entropias_Temperaturas(1,:,Pos);
Melhores_temperaturas = Matriz_Entropias_Temperaturas(2,:,Pos);

%Otimizar título

titulos_saida = 100*Matriz_Pressoes_Temp_reaq_Eficiencia_x14(:,end);
[Maior_titulo,Pos_2] = max(titulos_saida);
Melhor_config2 = Matriz_Pressoes_Temp_reaq_Eficiencia_x14(Pos_2,:)
            
%Curva de saturação

T = Sat.Temperature_C_;
s = Sat.Entropy_J_g_K_;

%Plot do diagrama T-s

hold on

plot(s,T,"--")
plot(Melhores_entropias,Melhores_temperaturas,".","Color","black")

title("Diagrama T-s")
xlabel("Entropia específica (kJ/kg.K)")
ylabel("Temperatura (°C)")

% Título na saída X Eficiência

[titulos_ord,ordem] = sort(titulos_saida);
Ef_ord = Eficiencias(ordem);

figure
plot(titulos_ord,Ef_ord,".")
title("Título da saída da TB X Eficiência do ciclo")
xlabel("Título (%)")
ylabel("Eficiência (%)")

toc