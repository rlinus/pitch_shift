classdef misc
%-------------------------------------------------------------------------
% Class name : Utils
%-------------------------------------------------------------------------
% Description : Contains misc utilities
%               * normalize
%               * normalizeMat
%-------------------------------------------------------------------------
%   Antoine Liutkus 2010
%-------------------------------------------------------------------------
    methods(Static)
        function V = normalize(V)
        % Normalizes V whether it be vector or matrix
            s = size(V);
            V = reshape(V(:)/sum(abs(V(:))),s);
        end
        
        function M = normalizeMat(M)
        % Normalizes each column of matrix M
            M = M - min(M(:));
            M = M./(ones(size(M,1),1)*sum(abs(M),1));
        end
        
        function V = nnrandn(varargin)
            V = abs(randn((varargin{:})));
        end
        
        function M = flip(M)
            M = flipud(fliplr(M));
        end
        
        function V =logistic(x)
            V = 1./(1+exp(-x));
        end
        
        function W = smoothWindow(N,Overlap)
            W = ones(N,1);
            if (Overlap < 0)||(Overlap>0.5)
                disp('SmoothWindow : warning, Overlap paremeter must be non negative and strictly < 0.5');
                Overlap = max(0,min(1,Overlap));
            end
            IOverlap = round(N*Overlap);
            W(1:IOverlap) = sin(linspace(0,pi/2,IOverlap)).^2;
            W(end-IOverlap+1:end) = sin(linspace(pi/2,0,IOverlap)).^2;
        end
        
        function d= divB(M1,M2,B)
            M1 = M1(:);
            M2 = M2(:);
            if nargin < 3
                B = 0; %default= itakura saito
            end
            switch B
                case 0
                    d = sum(M1./M2 - log(M1./M2) - 1);
                case 1
                    d = sum(M1.*(log(M1)-log(M2)) + (M1 - M2));
                otherwise
                    d = 1/B/(B-1)*sum(M1.^B + (B-1)*M2.^B - B*M1.*(M2.^(B-1))  );
            end
        end

        function NewMatrix = smooth(Matrix, Width)
            [N,T] = size(Matrix);
            FHorizon = round(2*Width);
            K = (-FHorizon:FHorizon);
            W = misc.normalize(exp(-K.^2/Width^2));
            LF = length(W);
            ZP = zeros(LF,T);
            NewMatrix = [ZP ; Matrix ; ZP];


            NewMatrix = filter(W, 1, NewMatrix,[],1);
            NewMatrix = NewMatrix(LF+FHorizon+(1:N),:);
        end        
        
        function Sc = completeSpectrum(S,nfft)
            nfftUtil = size(S,1);
            dimensions = size(S);
            Sc = zeros([nfft dimensions(2:end)]);
            dimensionString = repmat(',:',1,length(dimensions)-1);
            eval(['Sc(1:nfftUtil' dimensionString ') = S;']);
            eval(['Sc(end:-1:(end-nfftUtil+1)' dimensionString ') = conj(S);']);
        end
    end
end
    


