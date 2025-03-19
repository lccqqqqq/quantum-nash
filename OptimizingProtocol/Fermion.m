classdef Fermion
    
    properties
       ind; % label the sites
       dag; % distinguish c from cdaggers
       pref; % allow for scalar prefactors
    end
    
    methods
        function obj = Fermion(ind,dag,pref)
            obj.ind = ind;
            obj.dag = dag;
            obj.pref = pref;
            
            % check that the object is legal
            % ind must be a 1*n vector
            if and(isnumeric(ind),isnumeric(dag))
                if and(size(ind)==size(dag),or(size(ind,1)==1,size(ind,1)==0))
                    if numel(dag) == numel(dag(dag==0)) + numel(dag(dag==1))
                        return
                    else
                        error('dags must contain logical values')
                    end
                else
                    error('sizes does not match');
                end
            else
                error('[ind] and [dag] must be numeric');  
            end
            
        end
        
        
        
        function objs1 = cont(objs)
            % simplify individual terms
            objs1 = objs;
            
            for i = 1:numel(objs)
                f = objs(i);
                f1 = Fermion([],[],1);
%                 j = 1;
%                 while and(LengthReduction,j<=max(f.ind))
                if ~isempty(f.ind)
                     for j = unique(f.ind) % j:current fermion species
                        cc = []; % indices for creation ops dag=0
                        aa = []; % indices for annihilation ops dag=1
                        for l = 1:numel(f.ind)
                            if f.ind(l) == j
                                if f.dag(l) == 0
                                    aa = [aa,l];
                                elseif f.dag(l) == 1
                                    cc = [cc,l];
                                end
                            end       
                        end

                        if abs(numel(aa)-numel(cc)) > 1
                            f1.pref = 0;
                            break
                        elseif numel(aa)-numel(cc) == 1
                            indArray = [reshape([aa(1:end-1);cc],[1,2*numel(cc)]),aa(end)];
                            if issorted(indArray)
                                % follow fermion contraction rules
                                parity = sum(indArray)-sum(1:numel(indArray));
                                f.pref = f.pref*(-1)^parity;
                                f.dag(f.ind==j) = [];
                                f.ind(f.ind==j) = [];
                                f1.ind = [f1.ind,j];
                                f1.dag = [f1.dag,0];
                            else
                                f1.pref = 0;
                                break
                            end
                        elseif numel(cc)-numel(aa) == 1
                            indArray = [reshape([cc(1:end-1);aa],[1,2*numel(aa)]),cc(end)];
                            if issorted(indArray)
                                % follow fermion contraction rules
                                parity = sum(indArray)-sum(1:numel(indArray));
                                f.pref = f.pref*(-1)^parity;
                                f.dag(f.ind==j) = [];
                                f.ind(f.ind==j) = [];
                                f1.ind = [f1.ind,j];
                                f1.dag = [f1.dag,1];
                            else
                                f1.pref = 0;
                                break
                            end
                        else % elseif numel(cc) == numel(aa)
                            indArray1 = reshape([cc;aa],[1,2*numel(aa)]);
                            indArray2 = reshape([aa;cc],[1,2*numel(aa)]);
                            if issorted(indArray1) % creation on the left
                                parity = sum(indArray1)-sum(1:numel(indArray1));
                                f.pref = f.pref*(-1)^parity;
                                f.dag(f.ind==j) = [];
                                f.ind(f.ind==j) = [];
                                f1.ind = [f1.ind,j,j];
                                f1.dag = [f1.dag,1,0];
                            elseif issorted(indArray2)
                                parity = sum(indArray2)-sum(1:numel(indArray2)); 
                                f.pref = f.pref*(-1)^parity;
                                f.dag(f.ind==j) = [];
                                f.ind(f.ind==j) = [];
                                f1.ind = [f1.ind,j,j];
                                f1.dag = [f1.dag,0,1];
                            else
                                f1.pref = 0;
                                break
                            end
                        end

                        f1.pref = f.pref;
                     end
                else
                    % when the fermion is actually a scalar
                    f1.pref = f.pref;
                end
                
                if f1.pref == 0
                    objs1(i) = Fermion([],[],0);
                else
                    objs1(i) = f1;
                end
            end
            EmptyInds = [];
            for i = 1:numel(objs1)
                f = objs1(i);
                if f.pref == 0
                    EmptyInds = [EmptyInds,i];
                end
            end
            objs1(EmptyInds) = [];
        end
        
        function [tf,objs1] = isconted(objs)
            objs1 = cont(objs);
            tf = isequal(objs,objs1);
        end
        
        function objs1 = simplify(objs)
            objs1 = Fermion([],[],0); % empty fermion as initialization
            
            objs = cont(objs);
            tf = zeros(numel(objs));
            for i = 1:numel(objs)
                f = objs(i);
                for j = 1:numel(objs)
                    g = objs(j);
                    tf(i,j) = issimilar(f,g);
                end
            end
            
            for i = 1:numel(objs)
                if ~isequal(tf(:,i),zeros(numel(objs),1))
                    objs0 = objs(find(tf(:,i)));
                    pref0 = 0;
                    f0 = objs0(1);
                    for j = 1:numel(objs0)
                        f1 = objs0(j);
                        pref0 = pref0 + f1.pref;
                    end
                    % update tf
                    tf(find(tf(:,i)),:) = 0;
                    tf(:,find(tf(:,i))) = 0;
                    
                    % 
                    objs1 = [objs1,Fermion(f0.ind,f0.dag,pref0)];
                end
            end
            
            % check for zeros
            EmptyInds = [];
            for i = 1:numel(objs1) 
                f = objs1(i);
                if f.pref == 0
                    EmptyInds = [EmptyInds,i];
                end
            end
            objs1(EmptyInds) = [];
            
            if isempty(objs1)
                objs1 = Fermion([],[],0);
            end
            

            
        end
        
        function [objs2] = normal(objs)
            objs1 = Fermion([],[],0);
            for i = 1:numel(objs)
                f = objs(i);
                Nfset = Fermion([],[],1);
                
                if ~isempty(f.ind)
                    
                    [A,Aind] = unique(f.ind);
                    [~,ord] = sort(Aind);
                    
                    for j = A(ord)
                        ind1 = j*ones(1,sum(f.ind==j)); % can be 1 or 2
                        dag1 = f.dag(f.ind==j);
                        fs = Fermion(ind1,dag1,1);
                        if and(sum(fs.ind==j)>1,fs.dag(1)==0) % if config is c1*c1T
                            fs.dag = fliplr(fs.dag);
                            fs = [Fermion([],[],1),0-fs];
                        end 
                        Nfset = Nfset*fs;                      
                    end
                    Nfset = Nfset*f.pref;
                else
                    Nfset = f;
                end
                objs1 = [objs1,Nfset];
            end
            objs1(1) = [];
            
            objs2 = Fermion([],[],0);
            for i = 1:numel(objs1)

                f = objs1(i);
                if ~isempty(f.ind)
                    indices = 1:numel(f.ind);
                    % parity factor for normal ordering
                    p = 1;
                    p = p*(-1)^(sum(indices(f.dag==1))-sum(1:sum(f.dag==1)));
                    % parity factor for reverse order in the annihilation ops
                    p = p*(-1)^((sum(f.dag==0)-1)*sum(f.dag==0)/2);

                    cc = f.ind(f.dag==1);
                    aa = f.ind(f.dag==0);

                    ind1 = [cc,fliplr(aa)];
                    dag1 = [ones(size(cc)),zeros(size(aa))];
                    pref1 = f.pref*p;

                    objs2 = [objs2,Fermion(ind1,dag1,pref1)];
                else
                    objs2 = [objs2,objs1(i)];
                end
            end
            objs2(1) = [];
            
            objs2 = simplify(objs2);
        end
        
        
        

        
        function [obj] = plus(obj1,obj2)
            if and(isa(obj1,'Fermion'),isa(obj2,'Fermion'))
                obj = [obj1,obj2];
            elseif and(isa(obj1,'Fermion'),isa(obj2,'double'))
                obj = [obj1,Fermion([],[],obj2)];
            elseif and(isa(obj1,'double'),isa(obj2,'Fermion'))
                obj = [obj2,Fermion([],[],obj1)];
            end
            obj = normal(obj);
            
        end
        
        function [obj] = minus(obj1,obj2)
            if and(isa(obj1,'Fermion'),isa(obj2,'Fermion'))
                obj2m = obj2; % copy another Fermion chain with the same size
                for i = 1:numel(obj2)
                    f = obj2(i);
                    f.pref = -f.pref;
                    obj2m(i) = f;
                end
                obj = plus(obj1,obj2m);
            elseif and(isa(obj1,'Fermion'),isa(obj2,'double'))
                obj = plus(obj1,-obj2);
            elseif and(isa(obj1,'double'),isa(obj2,'Fermion'))
                obj2m = obj2; % copy another Fermion chain with the same size
                for i = 1:numel(obj2)
                    f = obj2(i);
                    f.pref = -f.pref;
                    obj2m(i) = f;
                end
                obj = plus(obj1,obj2m);
            end
            
            obj = simplify(obj);

        end
        
        function [obj] = mtimes(obj1,obj2)
            if and(isa(obj1,'Fermion'),isa(obj2,'Fermion'))
                obj = Fermion([],[],0);
                for i = 1:numel(obj1)
                    f = obj1(i);
                    for j = 1:numel(obj2)
                        g = obj2(j);
                        ind1 = [f.ind,g.ind];
                        dag1 = [f.dag,g.dag];
                        pref1 = f.pref*g.pref;

                        fg = Fermion(ind1,dag1,pref1);
                        obj = [obj,fg];
                    end
                end
            elseif and(isa(obj1,'Fermion'),isa(obj2,'double'))
                obj = obj1;
                for i = 1:numel(obj1)
                    f = obj1(i);
                    f.pref = f.pref*obj2;
                    obj(i) = f;
                end
            elseif and(isa(obj1,'double'),isa(obj2,'Fermion'))
                obj = obj2;
                for i = 1:numel(obj2)
                    f = obj2(i);
                    f.pref = f.pref*obj1;
                    obj(i) = f;
                end
            end
            obj = simplify(obj);
        end

        % need to check for identical elements
        function [tf] = isequal(objs1,objs2)
            % assume they are already simplified
            if numel(objs1)==numel(objs2)
                tf = true;
                for i = 1:numel(objs1)
                    f = objs1(i);
                    g = objs2(i);
                    tf = and(tf,isequal(f.ind,g.ind));
                    tf = and(tf,isequal(f.dag,g.dag));
                    tf = and(tf,isequal(f.pref,g.pref));
                end
            else
                tf = false;
            end
        end
        
        % check for similar terms
        function [tf] = issimilar(objs1,objs2)
            
            if numel(objs1)==numel(objs2)
                tf = true;
                for i = 1:numel(objs1)
                    f = objs1(i);
                    g = objs2(i);
                    tf = and(tf,isequal(f.ind,g.ind));
                    tf = and(tf,isequal(f.dag,g.dag));
                end
            else
                tf = false;
            end
           
        end
        
        function str = repr(obj)
            str = [];
            
            for i = 1:numel(obj)
                f = obj(i);
                
                
                
                if f.pref ~= 1
                    if f.pref < 0
                        str = [str,'(-',num2str(-f.pref),')'];
                    else
                        str = [str,num2str(f.pref)];
                    end
                end
                    
                
                if ~isempty(f.ind)
                    for j = 1:numel(f.ind)
                        str = [str,'c',num2str(f.ind(j))];
                        if f.dag(j) == 1
                            str = [str,'T'];
                        end
                        
                        if j < numel(f.ind)
                            str = [str,'*'];
                        end
                    end
                else
                    str = [str,num2str(f.pref)];
                end
                
                if i < numel(obj)
                    str = [str,' + '];
                end
            end
            
            
   
        end
        
        function obj = acomm(obj1,obj2)
            obj = obj1*obj2 + obj2*obj1;
            obj = normal(obj);
        end
            
        function obj = comm(obj1,obj2)
            obj = obj1*obj2 - obj2*obj1;
            obj = normal(obj);
        end
        
        function obj1 = hc(obj)
            obj1 = obj;
            for i = 1:numel(obj)
                f = obj(i);
                if ~isempty(f.ind)
                    fs = Fermion(fliplr(f.ind),1-fliplr(f.dag),conj(f.pref));
                    obj1(i) = fs;
                else
                    obj1(i) = Fermion([],[],conj(f.pref));
                end
            end
            obj1 = normal(obj1);
        end
        
        function m = vev(obj)
            obj = normal(obj);
            AllTermsNotEmpty = true;
            for i = 1:numel(obj)
                f = obj(i);
                AllTermsNotEmpty = and(AllTermsNotEmpty,~isempty(f.ind));
                if isempty(f.ind)
                    m = f.pref;
                    break   
                end
            end
            
            if AllTermsNotEmpty
                m = 0;
            end
            
        end
                    
    end
   
    methods (Static)
        % create multiple fermions
        function objs = plural(inds,dags,prefs)
            
            objs = Fermion([],[],0);
            for i = 1:numel(inds)
                f = Fermion(inds{i},dags{i},prefs(i));
                objs = [objs,f];
            end
            objs(1) = [];
            
        end
    end
end

