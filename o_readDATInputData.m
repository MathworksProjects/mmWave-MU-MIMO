function [ conf ] = o_readDATInputData( filePath, arguments )
%readDATInputDATA Read DAT Input Data File
%   Reads the DAT-formatted input data file and returns a struct with the
%   arguments and values read from the file, if found.
%   For arguments being arrays, the size may be indicated as the name of a
%   STRICTLY previously read variable: meaning that it has to appear before
%   in the arguments list
    lines = readLines(filePath);
    for i = 1:length(arguments)
        a = arguments{i};
        index = strncmp(lines,a.name,length(a.name));
        count = sum(index);
        if count == 1
            n = find(index == 1);
            l = lines{n};
            while ~contains(l, ';')
                l = [l,':',lines{n+1}]; %#ok<AGROW>
                n = n + 1;
            end
            if strcmp(a.type,'int')
                conf.(a.name) = readInt(l,a.name);
            elseif strcmp(a.type,'float')
                conf.(a.name) = readFloat(l,a.name);
            elseif strcmp(a.type,'char')
                conf.(a.name) = readChar(l,a.name);
            elseif strcmp(a.type,'bool')
                conf.(a.name) = readChar(l,a.name);
                if (strcmp(conf.(a.name),'true') || ...
                    strcmp(conf.(a.name),'True') || ...
                    strcmp(conf.(a.name),'T') || ...
                    strcmp(conf.(a.name),'t'))
                    conf.(a.name) = true;
                else
                    conf.(a.name) = false;
                end 
            elseif strcmp(a.type,'interval')
                conf.(a.name) = readInterval(l,a.name);
            elseif strcmp(a.type,'intArray')
                % a.size usually will be an integer, but it can also be a
                % string, meaning the name of a variable previously read
                if ischar(a.size)
                    conf.(a.name) = readArray(l,a.name,'int',conf.(a.size));
                elseif floor(a.size) == a.size
                    conf.(a.name) = readArray(l,a.name,'int',a.size);
                else
                    fprintf('The size of an array must be an integer\n');
                end
            elseif strcmp(a.type,'floatArray')
                if ischar(a.size)
                    conf.(a.name) = readArray(l,a.name,'float',conf.(a.size));
                elseif floor(a.size) == a.size
                    conf.(a.name) = readArray(l,a.name,'float',a.size);
                else
                    fprintf('The size of an array must be an integer\n');
                end
            elseif strcmp(a.type,'floatMatrix')
                if length(a.size) ~= 2
                    fprintf('The size of a matrix must be a 2-element array\n');
                    return
                end
                if ischar(a.size{1}) && ischar(a.size{2})
                    conf.(a.name) = readMatrix(l,a.name,'float',conf.(a.size{1}),conf.(a.size{2}));
                elseif floor(a.size{1}) == a.size{1} && ischar(a.size{2})
                    conf.(a.name) = readMatrix(l,a.name,'float',a.size{1},conf.(a.size{2}));
                elseif floor(a.size{2}) == a.size{2} && ischar(a.size{1})
                    conf.(a.name) = readMatrix(l,a.name,'float',conf.(a.size{1}),a.size{2});
                elseif floor(a.size{1}) == a.size{1} && floor(a.size{2}) == a.size{2}
                    conf.(a.name) = readMatrix(l,a.name,'float',a.size{1},a.size{2});
                else
                    fprintf('The column and row size of a matrix must be an integer\n');
                    return
                end
            elseif strcmp(a.type,'intMatrix')
                if length(a.size ~= 2) %#ok<ISMT>
                    fprintf('The size of a matrix must be a 2-element array\n');
                    return
                end
                if ischar(a.size{1}) && ischar(a.size{2})
                    conf.(a.name) = readMatrix(l,a.name,'int',conf.(a.size{1}),conf.(a.size{2}));
                elseif floor(a.size{1}) == a.size{1} && ischar(a.size{2})
                    conf.(a.name) = readMatrix(l,a.name,'int',a.size{1},conf.(a.size{2}));
                elseif floor(a.size{2}) == a.size{2} && ischar(a.size{1})
                    conf.(a.name) = readMatrix(l,a.name,'int',conf.(a.size{1}),a.size{2});
                elseif floor(a.size{1}) == a.size{1} && floor(a.size{2}) == a.size{2}
                    conf.(a.name) = readMatrix(l,a.name,'int',a.size{1},a.size{2});
                else
                    fprintf('The column and row size of a matrix must be an integer\n');
                    return
                end
            else
                fprintf('Unknown type %s\n',a.type)
                return
            end
            
        elseif count > 1
            fprintf('Argument %s duplicated\n',a.name)
            return
        else
            fprintf('Argument %s missing\n',a.name)
        end
        
        
    end
        

end

function [ arg ] = readInt( string, argument )
    nls = '\r\n';
    formatSpec = [argument,' = %d;',nls];
    arg = sscanf(string, formatSpec, 1);
end

function [ arg ] = readFloat( string, argument )
    nls = '\r\n';
    formatSpec = [argument,' = %f;',nls];
    arg = sscanf(string, formatSpec, 1);
end

function [ arg ] = readChar( string, argument )
    nls = '\r\n';
    formatSpec = [argument,' = %s;',nls];
    arg = sscanf(string, formatSpec, 1);
    arg(regexp(arg,'['';]'))=[];
end

function [ arg ] = readInterval( string, argument )
    nls = '\r\n';
    formatSpec = [argument,' = %s;',nls];
    arg = sscanf(string, formatSpec, 1);
    arg(regexp(arg,'['';]'))=[];
    arg = eval(arg); % The interval is expressed in MATLAB format
end

function [ arg ] = readArray( string, argument, type, size )
    nls = '\r\n';
    formatSpec = [argument,' = [ ',];
    if strcmp(type,'int')
        formatSpec = [formatSpec,repmat('%d ',1,size),'];',nls];
    elseif strcmp(type,'float')
        formatSpec = [formatSpec,repmat('%f ',1,size),'];',nls];
    end
    arg = sscanf(string, formatSpec, size)';
end

function [ arg ] = readMatrix( string, argument, type, Nx, Ny)
    arg = zeros(Nx,Ny);
    formatSpec = [argument,' = [:'];
    [~,~,~,NEXTINDEX] = sscanf(string, formatSpec, 1);
    lines = strsplit(string(NEXTINDEX:end),':');
    if strcmp(type,'int')
        formatSpec = ['[ ',repmat('%d\t',1,Ny),']'];
    elseif strcmp(type,'float')
        formatSpec = ['[ ',repmat('%f\t',1,Ny),']'];
    end
    for i = 1:Nx
        r = sscanf(lines{i}, formatSpec);
        if isempty(r)
            arg = [];
            return
        else
            arg(i,:) = r;
        end
    end
end

function [ lines ] = readLines( filePath )
    file = fopen(filePath,'r');
    sizeRead = 0;
    lines = cell(0);
    while sizeRead > -1
        r = fgets(file); % Transpose needed!
        if r == -1
            sizeRead = -1;
        else
            lines{end+1} = strtrim(r); %#ok<AGROW>
        end
    end
    fclose(file);
end

