function aa = cod2ami(seq)
if numel(seq) ~= 3
    error('Only 3 base pairs / string')
end

%% Turn all letters into capitals
seq = upper(seq);


switch seq(1)
    case 'T'
        switch seq(2)
            %% T T
            case'T' 
                switch seq(3)
                    case 'T' 
                        aa = 'phe';
                    case 'C'
                        aa = 'phe';
                    case 'A'
                        aa = 'leu';
                    case 'G'
                        aa = 'leu';
                    otherwise
                        error('You put garbage in, you get garbage out');
                end
                %% T C 
            case  'C'
                switch seq(3)
                    case 'T'
                        aa = 'ser';
                    case 'C'
                        aa = 'ser';
                    case 'A'
                        aa = 'ser';
                    case 'G'
                        aa = 'ser';
                    otherwise
                        error('You put garbage in, you get garbage out');
                end
                %% T A
            case 'A'
                switch seq(3)
                    case 'T'
                        aa = 'tyr';
                    case 'C'
                        aa = 'tyr';
                    case 'A'
                        aa = 'STOP';
                    case 'G'
                        aa = 'STOP';
                    otherwise
                        error('You put garbage in, you get garbage out');
                end
                %% T G
            case 'G'
                switch seq(3)
                    case 'T'
                        aa = 'cys';
                    case 'C'
                        aa = 'cys';
                    case 'A'
                        aa = 'STOP';
                    case 'G'
                        aa = 'trp';
                    otherwise
                        error('You put garbage in, you get garbage out');
                end
            otherwise
                error('You put garbage in, you get garbage out');
        end
    case 'C'
         switch seq(2)
            %% C T
            case 'T'
                switch seq(3)
                    case 'T'
                        aa = 'leu';
                    case 'C'
                        aa = 'leu';
                    case 'A'
                        aa = 'leu';
                    case 'G'
                        aa = 'leu';
                    otherwise
                        error('You put garbage in, you get garbage out');
                end
                %% C C 
            case 'C'
                switch seq(3)
                    case 'T'
                        aa = 'pro';
                    case 'C'
                        aa = 'pro';
                    case 'A'
                        aa = 'pro';
                    case 'G'
                        aa = 'pro';
                    otherwise
                        error('You put garbage in, you get garbage out');
                end
                %% C A
            case 'A'
                switch seq(3)
                    case 'T'
                        aa = 'his';
                    case 'C'
                        aa = 'his';
                    case 'A'
                        aa = 'gin';
                    case 'G'
                        aa = 'gin';
                    otherwise
                        error('You put garbage in, you get garbage out');
                end
                %% C G
            case 'G'
                switch seq(3)
                    case 'T'
                        aa = 'arg';
                    case 'C'
                        aa = 'arg';
                    case 'A'
                        aa = 'arg';
                    case 'G'
                        aa = 'arg';
                    otherwise
                        error('You put garbage in, you get garbage out');
                end
            otherwise
                error('You put garbage in, you get garbage out');
        end
    case 'A'
        switch seq(2)
            %% A T
            case 'T'
                switch seq(3)
                    case 'T'
                        aa = 'ile';
                    case 'C'
                        aa = 'ile';
                    case 'A'
                        aa = 'ile';
                    case 'G'
                        aa = 'met';
                    otherwise
                        error('You put garbage in, you get garbage out');
                end
                %% A C 
            case 'C'
                switch seq(3)
                    case 'T'
                        aa = 'thr';
                    case 'C'
                        aa = 'thr';
                    case 'A'
                        aa = 'thr';
                    case 'G'
                        aa = 'thr';
                    otherwise
                        error('You put garbage in, you get garbage out');
                end
                %% A A
            case 'A'
                switch seq(3)
                    case 'T'
                        aa = 'asn';
                    case 'C'
                        aa = 'asn';
                    case 'A'
                        aa = 'lys';
                    case 'G'
                        aa = 'lys';
                    otherwise
                        error('You put garbage in, you get garbage out');
                end
                %% A G
            case 'G'
                switch seq(3)
                    case 'T'
                        aa = 'ser';
                    case 'C'
                        aa = 'ser';
                    case 'A'
                        aa = 'arg';
                    case 'G'
                        aa = 'arg';
                    otherwise
                        error('You put garbage in, you get garbage out');
                end
            otherwise
                error('You put garbage in, you get garbage out');
        end
    case 'G'
        switch seq(2)
            %% G T
            case 'T'
                switch seq(3)
                    case 'T'
                        aa = 'val';
                    case 'C'
                        aa = 'val';
                    case 'A'
                        aa = 'val';
                    case  'G'
                        aa = 'val';
                    otherwise
                        error('You put garbage in, you get garbage out');
                end
                %% G C 
            case 'C'
                switch seq(3)
                    case 'T'
                        aa = 'ala';
                    case 'C'
                        aa = 'ala';
                    case 'A'
                        aa = 'ala';
                    case 'G'
                        aa = 'ala';
                    otherwise
                        error('You put garbage in, you get garbage out');
                end
                %% G A
            case 'A'
                switch seq(3)
                    case 'T'
                        aa = 'asp';
                    case 'C'
                        aa = 'asp';
                    case 'A'
                        aa = 'glu';
                    case 'G'
                        aa = 'glu';
                    otherwise
                        error('You put garbage in, you get garbage out');
                end
                %% G G
            case 'G'
                switch seq(3)
                    case 'T'
                        aa = 'gly';
                    case 'C'
                        aa = 'gly';
                    case 'A'
                        aa = 'gly';
                    case 'G'
                        aa = 'gly';
                    otherwise
                        error('You put garbage in, you get garbage out');
                end
            otherwise
                error('You put garbage in, you get garbage out');
        end
    otherwise
        error('You put garbage in, you get garbage out');
end