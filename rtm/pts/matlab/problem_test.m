sz=10;
A=2.*eye(sz,sz)+...
    [zeros(sz,1),[eye(sz-1,sz-1);zeros(1,sz-1)]]+...
    [[zeros(1,sz-1);eye(sz-1,sz-1)],zeros(sz,1)];
b=ones(sz,1);

A\b
