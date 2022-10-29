function [xsol, fsol, gsol] = knapsackProblem(W, c, w)
% 0-1背包问题，输入限重，各物品价值和重量的列向量，输出结果的序列x、总价值f和总质量g
    n = size(c, 1); % 根据c的行数确定物品数量，此题为20

    % 第一步：确定种群规模（只支持偶数），迭代次数，交叉概率，变异概率
    Np = 40; T = 100; pc = 0.8; pm = 0.3;
    % 第二步：随机生成初始种群，并且计算适应度函数f(x)
    X = rand(Np, n) > 0.5;      % 初始种群
    
    for t = 1 : T % 开始迭代
        [f, ~, J] = calculatef(W, X, c, w);     % 计算适应度函数
        % 第三步：利用轮盘赌法选择父代解，轮盘赌运行Np/2次
        parent = zeros(Np, n); offspring = parent;
        ff = f;     % 复制，用于修改
        for i = 1 : Np/2
            % 两轮轮盘赌，得到一对父代解
            j1 = RWS(ff); 
            parent(2*i-1, :) = X(j1, :);
            j2 = RWS(ff); 
            parent(2*i, :) = X(j2, :);
        
            % 第四步：随机数决定该对父代解是否进行交叉，单点交叉
            r = rand(1);               % 决定是否进行交叉操作
            if r < pc
                r = ceil(rand(1) * n); % 决定交叉点的随机数，ceil是向上取整
                % 交叉操作
                if r == 1
                    offspring(2*i-1, :) = parent(2*i, :);
                    offspring(2*i, :) = parent(2*i-1, :);
                else
                    offspring(2*i-1, :) = [parent(2*i-1, 1:r-1), parent(2*i, r:end)];
                    offspring(2*i, :) = [parent(2*i, 1:r-1), parent(2*i-1, r:end)];
                end
            else 
                offspring(2*i-1, :) = parent(2*i-1, :);
                offspring(2*i, :) = parent(2*i, :);
            end
        end
        
        % 第五步：变异
        mutation = rand(Np, n) < pm;    % Np*n的矩阵，由随机的0和1组成，为1的概率是pm
        offspring = abs(offspring - mutation);
        
        [f1, g1] = calculatef(W, X, c, w);
        [f2, g2] = calculatef(W, offspring, c, w);
        % 合并
        XX = [X;offspring];
        f = [f1; f2];
        g = [g1; g2];
        % 创造下一代种群
        ff = f;
        for i = 1 : Np
            % 储存最优的Np个解，作为下一代种群
            [~, j] = max(ff);   % 下一个最优解f值及其索引
            ff(j) = 0;          % 已遍历过的解设为0
            X(i, :) = XX(j, :);
        end
    end
    % 迭代结束，储存最终结果
    [f, g] = calculatef(W, X, c, w);
    [~, i] = max(f);
    xsol = X(i, :); fsol = f(i); gsol = g(i);
    % 不可行解的处理
    gg = 0;
    for j = J
        if (gg + xsol(j) * w(j) > W)  % 超重处理
            for jj = J(j:end)
                xsol(jj) = 0;
            end
            break                       %退出循环
        end
        gg = gg + xsol(j) * w(j) > W;
    end
end

function [f, g, J] = calculatef(W, X, c, w)
% 输入限重，种群，价值和重量，计算适应度函数f，以及总重g，同时返回c/w降序排列索引
    Np = size(X, 1);
    ratio = c./w;        % c/w比率
    % 不可行解的处理方法：解码方法
    ratio1 = ratio;             % 复制一份ratio用于循环中修改、使用
    J = [];                     % 用于索引储存
    for i = 1 : size(X, 2)
        % 该循环用于确定c/w比率降序排列后各元素的原索引
        [~, j] = max(ratio1);   % 下一个最大比率值及其索引
        ratio1(j) = 0;          % 已遍历过的值设为0
        J = [J, j];             % 按顺序储存索引
    end
    
    f = zeros(Np, 1); g = zeros(Np, 1);
    for k = 1 : Np
        % 该层循环用于处理第 k 个解，即第 k 行
        for j = J
            % 该层循环按顺序（J中已储存索引）选择物品，并得出g和f
            if (g(k) + X(k, j) * w(j) > W)  % 超重处理
                break                       % 退出循环
            end
            g(k) = g(k) + X(k, j) * w(j);   % 总重量
            f(k) = f(k) + X(k, j) * c(j);   % 总价值
        end
    end
end


function j = RWS(f)
% 轮盘赌选择，读取适应度函数列向量f，返回索引
    P = f / sum(f);             % 选择概率
    PP = [0; cumsum(P)];        % 概率的累积和
    r = rand(1);
    for j = 1 : size(f, 1)  
        if (r >= PP(j)) && (r <= PP(j+1))
            break
        end
    end
end