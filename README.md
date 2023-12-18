# Newton-Rapshon-method-non-linear-system
This is a Matlab implementation of Newton-Rapshon method that resolve an example of non linear system.

{█(e^(x+y)-1=0@e^(x-y)-1=0)┤

            %% HOMEWORK H1_3
            clc;
            clear; 
            close all;
            
            %% 1 - Varbs inizialization
            
            f1 = @(x,y) exp(x+y)-1; %f1
            
            f2 = @(x,y) exp(x-y)-1;%f2
            
            J = @(x,y) [exp(x+y) exp(x+y);
                        exp(x-y) -exp(x-y)]; %Jacobian
            
            x = input('Insert x0 = '); %values not so far from the true solution (0)
            y = input('Insert y0 = ');
            x0 = x;
            y0 = y;
            
            tol = 1e-4; 
            maxiter = 1e3;                                 
            h = [Inf;Inf];
             
            %varbs for the plots
            p = 0;
            hn = []; 
            hx = [];
            hy = [];
            Y = [];
            X = [];
            F1 = [];
            F2 = [];
            
            %% 2 - Newton-Raphson
            
            n = 0;
            while norm(h)>tol && n<=maxiter %cond exit
                n=n+1;
            
                if isnan(det(J(x,y))) || abs(det(J(x,y))) < tol %control convergence
                    p = 1;
                    disp('###############!WARNING!###############')
                    disp(['With input values x0 = ', num2str(x0),' e y0 = ', num2str(y0), ...
                          ' the algorithm not converges.']) 
                    disp(['x result = ' num2str(x)])        
                    disp(['y result = ' num2str(y)])      
                    disp(['Iterations = ' num2str(n)])         
                    disp('###############!WARNING!###############')
                    break                                      
                end                                           
                                            
                h = -J(x,y)\[f1(x,y);
                             f2(x,y)]; %increment
            
                hx(n) = h(1); %store varbs plots
                hy(n) = h(2);
                X(n) = x;
                Y(n) = y;
                F1(n) = f1(x,y);
                F2(n) = f2(x,y);
            
                x = x + h(1); %x_n+1
                y = y + h(2); %y_n+1
            end
            
            %% 3 - Plots
            
            %figure path convergence
            if p == 0
                figure(1) 
                hold on
                plot3(X,Y,F1,'-ko','LineWidth',1.5)
                grid on
                hold on
                plot3(X,Y,F2,'-rx','LineWidth',1.5)
                hold on
                plot(X(1),Y(1),'b','Marker','square', 'MarkerSize', 18); %starts
                hold on
                plot(X(end),Y(end),'g','Marker','square','MarkerSize', 18); %numerical solution
                hold off
                title(['Algorithm convergence (iterations = ', num2str(n) ')'],'FontSize',18);
                xlabel('x','FontSize',18)
                ylabel('y','FontSize',18)
                zlabel('numerical solution','FontSize',18)
                lgd2 = legend('Pathway of f_1(x,y) = e^x^+^y-1', ...
                       'Pathway of f_2(x,y) = e^x^-^y-1', ...
                       ['Start x_0 =' num2str(X(1)) ' e y_0 =' num2str(Y(1))], ...
                       ['Numerical solution  x =' num2str(x) ' y =' num2str(x)],'Location','northwest');
                set(lgd2, 'FontSize', 14)
                view(0, 90);
            
            %increment pathway
                figure(2) 
                plot3(X,Y,hx,'-ko','LineWidth',1.5)
                grid on
                hold on
                plot3(X,Y,hy,'-rx','LineWidth',1.5)
                hold on
                plot3(X(1),Y(1),hx(1),'Marker','square', 'MarkerSize', 12,'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g') %punto iniziale
                hold on                                                
                plot3(X(1),Y(1),hy(1),'Marker','square', 'MarkerSize', 12,'MarkerEdgeColor','r',...
                    'MarkerFaceColor','g') %punto iniziale
                title(['Increment of x and y (iterations = ', num2str(n) ')'], ...
                       'FontSize',18);
                xlabel('x','FontSize',18);
                ylabel('y','FontSize',18);
                zlabel('h','FontSize',18)
                lgd3 = legend('x increment (h_x_n)','y increment (h_y_n)','start','Location','northwest');
                set(lgd3, 'FontSize', 14)
                view(0, 90);
            
            %disp risultati
                disp(' ')
                disp('###############-RESULTS-###############')
                disp(['Numerical solution x = ' num2str(x)]) %disp res
                disp(['Numerical solution y = ' num2str(y)])
                disp(['Iterations = ' num2str(n)])
            end
