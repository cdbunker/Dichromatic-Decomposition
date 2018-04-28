A = im2double(imread('5.png'));
A = imresize(A, 1);

ms = 0.6*ones(size(A,1), size(A,2));
md = 0.7*ones(size(A,1), size(A,2));
cs = [mean(mean(A(:,:,1))); mean(mean(A(:,:,2))); mean(mean(A(:,:,3)))];
k = 100;

[L,K] = superpixels(A,k);

figure
imshow(A)
title('Original Image');

figure
BW = boundarymask(L);
imshow(imoverlay(A,BW,'cyan'),'InitialMagnification',67)
title('Image showing superpixel boundaries');

cdDense = zeros(size(A),'like',A);
idx = label2idx(L);
numRows = size(A,1);
numCols = size(A,2);
for labelVal = 1:K
    redIdx = idx{labelVal};
    greenIdx = idx{labelVal}+numRows*numCols;
    blueIdx = idx{labelVal}+2*numRows*numCols;
    cdDense(redIdx) = mean(A(redIdx));
    cdDense(greenIdx) = mean(A(greenIdx));
    cdDense(blueIdx) = mean(A(blueIdx));
    cd(1,labelVal) = mean(A(redIdx));
    cd(2,labelVal) = mean(A(greenIdx));
    cd(3,labelVal) = mean(A(blueIdx));
    mat{labelVal} = [redIdx, md(redIdx), A(redIdx), A(greenIdx), A(blueIdx), ms(redIdx)];
end

figure
imshow(cdDense,'InitialMagnification',67)
title('Average pixel value within each superpixel');

for count = 1:10
    %tic
    fprintf('Estimating Cd\n')
    for labelVal = 1:K
        s = 0;
        sMd = 0;
        info = mat{labelVal};
        %info is index, md, R, G, B, ms
        l = length(info);
        for i = 1:l
            idx = info(i,1);
            innerTerm = ms(idx)*cs - info(i, 3:5)';
            s = s + (md(idx)*innerTerm);
            sMd = sMd + md(idx).^2;
        end
        cd(:,labelVal) = s ./ sMd;
    end
    cd=abs(cd);
    cd=normc(cd);
    %toc
    
    fprintf('\nEstimating md\n')
    %tic
    for labelVal = 1:K
        info = mat{labelVal};
        %info is index, md, R, G, B, ms
        l = length(info);
        for i = 1:l
            idx = info(i,1);
            top = cd(:,labelVal)'*(ms(idx)*cs - info(i, 3:5)');
            bottom = cd(:,labelVal)'*cd(:,labelVal);
            md(idx) = top/bottom;
        end
    end
    md=abs(md);
    %toc
    
    fprintf('\nEstimating Cs\n')
    %tic
    s=0;
    for labelVal = 1:K
        info = mat{labelVal};
        %info is index, md, R, G, B, ms
        l = length(info);
        for i = 1:l
            idx = info(i,1);
            s = s + ms(idx)*(md(idx)*cd(:,labelVal) - info(i, 3:5)');
        end
        cs = s/sum(sum(ms.^2));
    end
    cs=abs(cs);
    cs=cs/norm(cs);
    
    %toc
    
    fprintf('\nEstimating ms\n')
    %tic
    for labelVal = 1:K
        info = mat{labelVal};
        %info is index, md, R, G, B, ms
        l = length(info);
        for i = 1:l
            idx = info(i,1);
            top = cs'*(md(idx)*cd(:,labelVal) - info(i, 3:5)');
            bottom = cs'*cs;
            ms(idx) = top/bottom;
        end
    end
    
    ms=abs(ms);
    ms=ms/max(ms(:));
    %toc
    
    cdDense2 = zeros(size(A));
    for i=1:(size(A,1)*size(A,2))
        group = L(i);
        [I,J] = ind2sub(size(A), i);
        cdDense2(I,J,:) = cd(:,group);
    end
    
    csDense2 = repmat(reshape(cs,1,1,3),size(A, 1), size(A,2));
    
    reconstructed = md.*cdDense2 + ms.*csDense2;
    error = mean(sum(sum((A-reconstructed).^2))/(size(A,1)*size(A,2)));
    fprintf('%1.10f\n', error)
end

cdDense2 = zeros(size(A));
for i=1:(size(A,1)*size(A,2))
    group = L(i);
    [I,J] = ind2sub(size(A), i);
    cdDense2(I,J,:) = cd(:,group);
end

csDense2 = repmat(reshape(cs,1,1,3),size(A, 1), size(A,2));

reconstructed = md.*cdDense2 + ms.*csDense2;
error = mean(sum(sum((A-reconstructed).^2))/(size(A,1)*size(A,2)));

cdDense2 = zeros(size(A));
for i=1:(size(A,1)*size(A,2))
    group = L(i);
    [I,J] = ind2sub(size(A), i);
    cdDense2(I,J,:) = cd(:,group);
end

csDense2 = repmat(reshape(cs,1,1,3),size(A, 1), size(A,2));

reconstructed = md.*cdDense2 + ms.*csDense2;
error = mean(sum(sum((A-reconstructed).^2))/(size(A,1)*size(A,2)))

figure
imshow(csDense2)
title('Cs');
figure
imshow(ms/(max(ms(:))))
title('ms');
figure
imshow(cdDense2)
title('Cd');
figure
imshow(md/(max(ms(:))))
title('md');

figure
imshow(reconstructed)
title('md * Cd + ms * Cs');
figure
imshow(md.*cdDense2)
title('md * Cd');
figure
imshow(ms.*csDense2)
title('ms * Cs');
