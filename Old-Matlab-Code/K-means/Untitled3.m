% tic tac toe

close all; clc; clear variables
while true
    board = magic(3)
    
    choose1 = randi(2);
    
    
    
    if choose1 == 1;
        pick = input('User goes first and is 0, choice: ');
        flag = 0;
    else
        pick = randi(9);
        disp(['Comp when first and is -1 picked: ',num2str(pick)]);
        flag =1;
    end
    
    [row, col] = find(board == pick);
    breakout = 0;
    if choose1 == 1;
        board(row,col) = 0
    else
        board(row,col) = -1
    end
    block = 0;
    winpick = 0;
    while true
        if flag == 1;
            pick = input('User choice: ');
            
        else
            if block == 1
                pick = newpick;
                if winpick ~= 0
                    pick = winpick;
                    winpick = 0;
                end
                block = 0;
            elseif pick == 5 && block ~= 1
                goco = randi(4);
                switch goco
                    case 1
                        pick = board(1,1);
                    case 2
                        pick = board(1,3);
                    case 3
                        pick = board(3,1);
                    case 4
                        pick = board(3,3);
                    otherwise
                end
            else
                pick = randi(9);
            end
        end
        [row, col] = find(board == pick);
        if isempty(row) || pick == 0 || pick == -1
            if flag == 1;
                disp('User choice is invalid, try picking something less stupid');
            end
        else
            clc
            if flag == 1
                board(row,col) = 0
                flag = 0;
            else
                disp(['Comp picked: ',num2str(pick)]);
                board(row,col) = -1
                flag = 1;
            end
        end
        
        for i = 1:3
            if isequal(board(i,:),[0,0,0]) || isequal(board(:,i),[0,0,0].')
                disp('User wins!');
                breakout = 1;
            end
            if isequal(board(i,:),[-1,-1,-1]) || isequal(board(:,i),[-1,-1,-1].')
                disp('Comp wins user sucks!');
                breakout = 1;
            end
        end
        if isequal(diag(board), [0,0,0].') || isequal(diag(fliplr(board)), [0,0,0].')
            disp('User wins!');
            breakout = 1;
        end
        if isequal(diag(board), [-1,-1,-1].') || isequal(diag(fliplr(board)), [-1,-1,-1].')
            disp('Comp wins!');
            breakout = 1;
        end
        [row,col] = find(board>0);
        if isempty(row)
            disp('Draw!');
            breakout = 1;
        end
        if breakout == 1
            break
        end
        choi = 0;
        if flag == 0
            [row, col] = find(board > 0); 
            bl_board = board;
            for j = 1:numel(row)
                bl_board(row(j), col(j)) = 0;
                for i = 1:3
                    if isequal(bl_board(i,:),[0,0,0]) || isequal(bl_board(:,i),[0,0,0].') || isequal(diag(bl_board), [0,0,0].') || isequal(diag(fliplr(bl_board)), [0,0,0].')
                        block = 1;
                        newpick = board(row(j),col(j));
                        choi = 1;
                    end
                end
                
                bl_board(row(j), col(j)) = -1;
                for i = 1:3
                    if isequal(bl_board(i,:),[-1,-1,-1]) || isequal(bl_board(:,i),[-1,-1,-1].') || isequal(diag(bl_board), [-1,-1,-1].') || isequal(diag(fliplr(bl_board)), [-1,-1,-1].')
                        block = 1;
                        winpick = board(row(j),col(j));
                        choi = 1;
                    end
                end
                bl_board = board;
%                 if choi == 1;
%                     break
%                 end
            end
        end
    end
    an = input('Play again? ' ,'s');
    if 1-strcmp(an,'y') && 1 - strcmp(an,'Y')
        break
    end
    clc
end