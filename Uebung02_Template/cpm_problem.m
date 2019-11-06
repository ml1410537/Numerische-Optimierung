function [A, b, c] = cpm_problem(activities, durations, predecessors)
%CPM_PROBLEM Create linear program for critical path method
%
% Vorlesung 'Numerische Optimierung und modellprädiktive Regelung'
%
% Template für Methode des kritischen Pfades
%
% Ergänzen Sie die '???'-Einträge
%
% ?bergabeargumente:
%
% activities    Aufgaben
% durations     Dauer der Aufgaben
% predecessors  Paare mit Vorgänger, Nachfolger

    % number of activities
    k = length(activities);
    % number of equality constraints
    p = k + size(predecessors, 1);
    % number of variables (s_i, T, slack)
    n = k + 1 + p;
    
    % cost vector
    c = zeros(n, 1);
    % to do
    c(k+1) = 1;
    % equality constraints
    A = zeros(p, n);
    b = zeros(p, 1);
    for i = 1 : k
        % to do
        A(i,i) = -1;
        A(i,k+1) = 1;
        A(i,k+1+i) = -1;
        b(i) = durations(i);
    end
    for idx = 1 : size(predecessors, 1)
       i = predecessors(idx, 1);
       j = predecessors(idx, 2);
       % to do
       A(k+idx,i) = -1;
       A(k+idx,j) = 1;
       A(k+idx,k+1+k+idx) = -1;
       b(k+idx) = durations(i);
    end
end
