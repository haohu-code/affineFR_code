function satProblem = loadSATLIB(filename)
    % Open the DIMACS CNF file
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file: %s', filename);
    end
    cleanupObj = onCleanup(@() fclose(fid));

    satProblem.clauses = {};
    satProblem.numVars = 0;
    satProblem.numClauses = 0;
    currentClause = [];
    foundProblemLine = false;
    lineNumber = 0;

    while ~feof(fid)
        line = fgetl(fid);
        lineNumber = lineNumber + 1;
        if ~ischar(line)
            break;
        end
        line = strtrim(line);
        
        % Skip empty lines or comments
        if isempty(line) || startsWith(line, 'c')
            continue;
        end
        
        % Parse the problem line: p cnf [vars] [clauses]
        if startsWith(line, 'p')
            parts = strsplit(line);
            if length(parts) < 4 || ~strcmp(parts{2}, 'cnf')
                error('Invalid DIMACS problem line in %s at line %d.', ...
                    filename, lineNumber);
            end
            satProblem.numVars = str2double(parts{3});
            satProblem.numClauses = str2double(parts{4});
            if ~isfinite(satProblem.numVars) || ...
                    ~isfinite(satProblem.numClauses) || ...
                    satProblem.numVars < 0 || satProblem.numClauses < 0 || ...
                    satProblem.numVars ~= floor(satProblem.numVars) || ...
                    satProblem.numClauses ~= floor(satProblem.numClauses)
                error('Invalid variable or clause count in %s at line %d.', ...
                    filename, lineNumber);
            end
            foundProblemLine = true;
            continue;
        end
        
        if ~foundProblemLine
            error('Clause data precedes the DIMACS problem line in %s at line %d.', ...
                filename, lineNumber);
        end

        % DIMACS CNF is a token stream: a zero terminates a clause, while
        % line breaks have no mathematical meaning. Thus, a clause may span
        % several lines and a single line may contain several clauses.
        nums = sscanf(line, '%d')';
        for k = 1:length(nums)
            literal = nums(k);
            if literal == 0
                satProblem.clauses{end+1} = currentClause;
                currentClause = [];
            else
                if abs(literal) > satProblem.numVars
                    error(['Literal %d exceeds the declared number of ' ...
                        'variables in %s at line %d.'], ...
                        literal, filename, lineNumber);
                end
                currentClause(end+1) = literal;
            end
        end
    end

    if ~foundProblemLine
        error('Missing DIMACS problem line in %s.', filename);
    end
    if ~isempty(currentClause)
        error('Unterminated clause at the end of %s.', filename);
    end
    if length(satProblem.clauses) ~= satProblem.numClauses
        error(['Clause-count mismatch in %s: header declares %d clauses, ' ...
            'but %d zero-terminated clauses were parsed.'], ...
            filename, satProblem.numClauses, length(satProblem.clauses));
    end
    
    % Optional: Convert to a numeric matrix if all clauses have same length (e.g., 3-SAT)
    if all(cellfun(@length, satProblem.clauses) == 3)
        satProblem.clauseMatrix = cell2mat(satProblem.clauses');
    end
end
