% Resistor Circuit Solver using Modified Nodal Analysis (MNA)
% Solves DC circuits with resistors and independent voltage/current sources.

function circuitSolverMNA()

  while true % Loop for re-running the program
      clearvars -except circuitSolverMNA; % Clear variables except the function itself
      clc; % Clear command window

      fprintf('--- 基于MNA的电阻电路求解器 ---\n');
      fprintf('请使用 NaN 表示未知值 "x"。\n');
      fprintf('电压方向定义为：行编号节点 -> 列编号节点 (V_col - V_row)。\n');
      fprintf('电流方向定义为：行编号节点 -> 列编号节点。\n\n');

      % --- 1. 获取用户输入 ---
      n = input('请输入节点个数 (n): ');
      if isempty(n) || ~isnumeric(n) || n < 2 || floor(n) ~= n
          fprintf('错误：节点个数必须是大于等于2的整数。\n');
          continue; % Restart the loop
      end

      fprintf('\n--- 输入节点连接矩阵 (n x n) ---\n');
      fprintf('仅需输入上三角部分 (i < j)，1表示连接，0表示断开。\n');
      fprintf('其余部分将自动处理。\n');
      Connectivity = getUserMatrixInput(n, 'Connectivity', [0, 1]);
      if isempty(Connectivity), continue; end % Restart if input failed

      fprintf('\n--- 输入节点间电阻矩阵 (n x n, 单位: Ohm) ---\n');
      fprintf('仅需输入上三角部分 (i < j)。输入电阻值，或用 NaN 表示未知/非电阻元件。\n');
      Resistance = getUserMatrixInput(n, 'Resistance', [0, Inf]); % Allow 0 and Inf, NaN for unknown/other
      if isempty(Resistance), continue; end

      fprintf('\n--- 输入节点间电压源矩阵 (n x n, 单位: Volt) ---\n');
      fprintf('仅需输入上三角部分 (i < j)。输入电压源值 (Vj - Vi)，或用 NaN 表示未知/无电压源。\n');
      Voltage = getUserMatrixInput(n, 'Voltage', [-Inf, Inf]); % Allow any voltage, NaN for unknown/none
      if isempty(Voltage), continue; end

      fprintf('\n--- 输入节点间电流源矩阵 (n x n, 单位: Ampere) ---\n');
      fprintf('仅需输入上三角部分 (i < j)。输入电流源值 (i -> j)，或用 NaN 表示未知/无电流源。\n');
      Current = getUserMatrixInput(n, 'Current', [-Inf, Inf]); % Allow any current, NaN for unknown/none
      if isempty(Current), continue; end

      % --- 2. 验证输入一致性 (基本检查) ---
      % Check if a branch is defined multiple ways (e.g., R and Vs both specified)
      inconsistent_branch = false;
      for r = 1:n
          for c = r+1:n
              if Connectivity(r, c) == 1
                  num_definitions = ~isnan(Resistance(r,c)) + ~isnan(Voltage(r,c)) + ~isnan(Current(r,c));
                  % Allow Resistance=0 (short) or Resistance=Inf (open) together with NaN for V/I
                  is_short_or_open = ~isnan(Resistance(r,c)) && (Resistance(r,c) == 0 || isinf(Resistance(r,c)));
                  if is_short_or_open && num_definitions > 1 && (~isnan(Voltage(r,c)) || ~isnan(Current(r,c)))
                      fprintf('警告: 节点 %d 和 %d 之间定义了短路/开路电阻，同时指定了电压/电流源。将忽略电阻值。\n', r, c);
                       Resistance(r,c) = NaN; % Prioritize V/I source over R=0 or R=Inf if specified
                       num_definitions = ~isnan(Resistance(r,c)) + ~isnan(Voltage(r,c)) + ~isnan(Current(r,c));
                  elseif ~is_short_or_open && num_definitions > 1
                      fprintf('错误: 节点 %d 和 %d 之间定义了多个元件 (电阻/电压源/电流源)。请只指定一个。\n', r, c);
                      inconsistent_branch = true;
                      break; % Exit inner loop
                  end
              elseif Connectivity(r, c) == 0 && (~isnan(Resistance(r,c)) || ~isnan(Voltage(r,c)) || ~isnan(Current(r,c)))
                   fprintf('警告: 节点 %d 和 %d 标记为未连接，但输入了元件值。将忽略这些值。\n', r, c);
                   Resistance(r,c) = NaN;
                   Voltage(r,c) = NaN;
                   Current(r,c) = NaN;
              end
          end
          if inconsistent_branch, break; end % Exit outer loop
      end
      if inconsistent_branch
           fprintf('条件不满足，请重新输入。\n\n');
           continue; % Restart main loop
      end


      % --- 3. 构建改进节点分析 (MNA) 矩阵 ---
      % 参考节点选择最后一个节点 n
      num_nodes = n - 1; % Number of non-reference nodes
      num_vs = 0; % Counter for voltage sources
      vs_indices = []; % Store info about voltage sources [row, col, index]

      % Identify voltage sources to determine matrix size
      for r = 1:n
          for c = r+1:n
              if Connectivity(r, c) == 1 && ~isnan(Voltage(r, c))
                  num_vs = num_vs + 1;
                  vs_indices(end+1, :) = [r, c, num_vs]; % Store node pair and VS index
              end
          end
      end

      matrix_size = num_nodes + num_vs;
      A = zeros(matrix_size, matrix_size); % MNA matrix
      b = zeros(matrix_size, 1);       % Right-hand side vector

      % Fill G (conductance) part and Is part (current sources)
      for r = 1:num_nodes % Iterate through non-reference nodes for KCL
          % Diagonal elements of G
          for c = 1:n % Check connection to all other nodes (including reference)
              node1 = r;
              node2 = c;
              % Ensure lower index comes first for accessing matrices
              if node1 > node2
                  tmp = node1;
                  node1 = node2;
                  node2 = tmp;
              end

              if node1 == node2 continue; end % Skip self

               if Connectivity(node1, node2) == 1
                   % Contribution from Resistors
                   if ~isnan(Resistance(node1, node2))
                       res_val = Resistance(node1, node2);
                       if res_val > 0 && ~isinf(res_val) % Standard resistor
                           conductance = 1 / res_val;
                           A(r, r) = A(r, r) + conductance; % Add to diagonal
                           if c <= num_nodes % If connected to another non-ref node
                               A(r, c) = A(r, c) - conductance; % Subtract from off-diagonal
                           end
                       elseif res_val == 0 % Short circuit (ideal wire) - Handle carefully
                           % MNA typically handles shorts via voltage source V=0
                           % Or by merging nodes - simpler here is to treat as Vs=0
                           % Let's add it as a V=0 source if not already a V source
                           is_vs = ~isnan(Voltage(node1, node2));
                           if ~is_vs
                               num_vs = num_vs + 1;
                               vs_indices(end+1, :) = [node1, node2, num_vs];
                               Voltage(node1, node2) = 0; % Treat short as V=0 source
                               % Resize matrices if needed (should pre-calculate size better)
                               if size(A,1) < num_nodes + num_vs
                                   new_size = num_nodes + num_vs;
                                   A(new_size, new_size) = 0;
                                   b(new_size, 1) = 0;
                               end
                               fprintf('信息: 节点 %d 和 %d 间电阻为0，视为零电压源处理。\n', node1, node2);
                           end
                            Resistance(node1,node2)=NaN; % Mark resistance as not primary element anymore
                       elseif isinf(res_val)
                           % Open circuit, conductance is 0, do nothing
                       end
                   end

                   % Contribution from Current Sources to b vector
                   if ~isnan(Current(node1, node2))
                       current_val = Current(node1, node2);
                       % Current flows node1 -> node2
                       if r == node1 % Current leaving node r
                           b(r) = b(r) - current_val;
                       elseif r == node2 % Current entering node r
                           b(r) = b(r) + current_val;
                       end
                   end
               end
          end
      end

      % Fill B, C, D parts (Voltage source contributions)
      current_vs_index = 0;
      for i = 1:size(vs_indices, 1)
          r = vs_indices(i, 1); % Node row index
          c = vs_indices(i, 2); % Node col index
          vs_num = vs_indices(i, 3); % Index of this voltage source current in X vector

          vs_val = Voltage(r, c);
          col_index = num_nodes + vs_num; % Column in A for this VS current

          % KCL equations (B matrix part and C^T part)
          if r <= num_nodes % If node r is not reference
              A(r, col_index) = A(r, col_index) + 1; % Current I_vs leaves node r
              A(col_index, r) = A(col_index, r) + 1; % For V_r in VS equation
          end
          if c <= num_nodes % If node c is not reference
              A(c, col_index) = A(c, col_index) - 1; % Current I_vs enters node c
              A(col_index, c) = A(col_index, c) - 1; % For V_c in VS equation
          end

          % Voltage source constraint equation (row = num_nodes + vs_num)
          b(num_nodes + vs_num) = vs_val; % Vc - Vr = Vs_val

          % D matrix part is usually zero for independent sources
          % A(num_nodes + vs_num, num_nodes + vs_num) = 0;
      end

      % --- 4. 求解线性方程组 AX = b ---
      % Check matrix condition before solving
      if rank(A) < matrix_size || cond(A) > 1e10 % Check rank and condition number
          fprintf('\n错误：电路条件不足或存在矛盾 (矩阵奇异或病态)。\n');
          fprintf('可能原因：\n');
          fprintf('- 电路未完全约束 (例如，浮空部分)。\n');
          fprintf('- 电压源环路或电流源节点不满足 KVL/KCL。\n');
          fprintf('- 输入错误。\n');
          fprintf('条件不满足，请重新输入。\n\n');
          continue; % Restart main loop
      end

      try
          X = A \ b; % Solve the system
      catch ME
          fprintf('\n错误：求解线性方程组失败。\n');
          fprintf('MATLAB 错误信息: %s\n', ME.message);
          fprintf('条件不满足，请重新输入。\n\n');
          continue; % Restart main loop
      end

      % --- 5. 提取并整理结果 ---
      NodeVoltages = zeros(n, 1);
      NodeVoltages(1:num_nodes) = X(1:num_nodes); % Assign solved voltages (V1 to Vn-1)
      NodeVoltages(n) = 0; % Reference node voltage

      VsCurrents = zeros(num_vs, 1);
      if num_vs > 0
           VsCurrents = X(num_nodes+1 : end); % Assign solved voltage source currents
      end

      % Initialize output matrices
      VoltOut = NaN(n, n); % Voltage difference V_j - V_i
      CurrOut = NaN(n, n); % Current flow i -> j
      ResOut = NaN(n, n);  % Equivalent Resistance V_ij / I_ij

      vs_current_map = containers.Map('KeyType','char','ValueType','double');
      for i = 1:size(vs_indices,1)
           key = sprintf('%d-%d', vs_indices(i,1), vs_indices(i,2));
           vs_current_map(key) = VsCurrents(i);
      end


      for r = 1:n
          for c = r+1:n
              if Connectivity(r, c) == 1
                  v_r = NodeVoltages(r);
                  v_c = NodeVoltages(c);
                  v_diff = v_c - v_r; % Voltage Vc - Vr
                  VoltOut(r, c) = v_diff;

                  current_ij = NaN; % Current from r to c

                  % Determine current based on original element type
                  if ~isnan(Resistance(r, c))
                      res_val = Resistance(r,c);
                       if res_val > 0 && ~isinf(res_val)
                            current_ij = v_diff / res_val;
                       elseif res_val == 0 % Short circuit (handled as Vs=0)
                           key = sprintf('%d-%d', r, c);
                           if isKey(vs_current_map, key)
                               current_ij = -vs_current_map(key); % Current r->c for Vs=Vc-Vr=0 source
                                                                 % MNA solves for current *leaving* r via VS branch
                           else
                               fprintf('警告：无法找到电阻为0的支路 (%d-%d) 的电流。\n', r, c);
                           end
                       elseif isinf(res_val)
                           current_ij = 0; % Open circuit
                       end
                       ResOut(r, c) = res_val; % Store original resistance

                  elseif ~isnan(Voltage(r, c)) % Voltage Source
                      key = sprintf('%d-%d', r, c);
                      if isKey(vs_current_map, key)
                           % MNA current I_Vs flows out of '+' node (c) into source, or out of node r into branch
                           % Our definition: CurrOut(r,c) is current r->c
                           % MNA variable I_Vs is current flowing through VS branch, positive direction needs care based on A matrix construction.
                           % In our A matrix: A(r, col_idx)=+1 (Ivs leaves r), A(c, col_idx)=-1 (Ivs enters c)
                           % So the solved Ivs is the current flowing r->c through the source.
                           current_ij = vs_current_map(key);
                      else
                           fprintf('警告：无法找到电压源 (%d-%d) 的电流。\n', r, c);
                      end
                      % Calculate equivalent resistance if current is non-zero
                      if ~isnan(current_ij) && abs(current_ij) > 1e-9 % Avoid division by zero
                         ResOut(r, c) = v_diff / current_ij;
                      elseif abs(v_diff) < 1e-9 && abs(current_ij) < 1e-9
                          ResOut(r, c) = NaN; % Indeterminate 0/0
                      else
                          ResOut(r, c) = Inf; % Essentially infinite resistance if current is zero
                      end

                  elseif ~isnan(Current(r, c)) % Current Source
                      current_ij = Current(r, c); % Current r -> c is fixed by source
                       % Calculate equivalent resistance
                       if abs(current_ij) > 1e-9
                           ResOut(r, c) = v_diff / current_ij;
                       elseif abs(v_diff) < 1e-9 && abs(current_ij) < 1e-9
                           ResOut(r, c) = NaN; % Indeterminate 0/0
                       else
                           ResOut(r, c) = Inf; % Infinite resistance if V!=0 but I=0 (unlikely for ideal source)
                                              % Or -Inf if V!=0 and I=0? Let's stick to Inf.
                       end
                  else
                     % Branch exists but wasn't defined R, Vs, or Is? Should not happen based on input checks
                     % Could happen if Connectivity=1 but R=NaN, V=NaN, I=NaN. Treat as open circuit.
                     current_ij = 0;
                     ResOut(r,c) = Inf;
                  end

                  CurrOut(r, c) = current_ij;

              end % if connected
          end % for c
      end % for r

      % --- 6. 显示结果 ---
      fprintf('\n--- 求解结果 ---\n');

      fprintf('\n节点电压 (相对节点 %d):\n', n);
      for i = 1:n
          fprintf('  V%d = %.4f V\n', i, NodeVoltages(i));
      end

      fprintf('\n节点间电压矩阵 (Vj - Vi, 上三角, 单位: V):\n');
      displayMatrixWithX(VoltOut);

      fprintf('\n节点间电流矩阵 (i -> j, 上三角, 单位: A):\n');
      displayMatrixWithX(CurrOut);

      fprintf('\n节点间等效电阻矩阵 (上三角, 单位: Ohm):\n');
      displayMatrixWithX(ResOut);


      % --- 7. 询问是否再次运行 ---
      fprintf('\n');
      rerun = input('是否要分析新的电路? (y/n): ', 's');
      if lower(rerun) ~= 'y'
          break; % Exit the while loop
      end
  end % End while true loop

  fprintf('程序结束。\n');

end % End main function

% --- 辅助函数：获取用户矩阵输入 ---
function matrix = getUserMatrixInput(n, matrixName, allowedRange)
  matrix = NaN(n, n); % Initialize with NaN
  validInput = false;
  while ~validInput
      fprintf('请输入 %s 矩阵的上三角部分 (行号 < 列号)。\n', matrixName);
      fprintf('用空格分隔列，用分号分隔行（如果需要多行输入）。\n');
      fprintf('对于 %d x %d 矩阵，您需要输入 %d 行（或更少，如果后面行为空）。\n', n, n, n-1);
      fprintf('示例 (3x3): val12 val13; val23\n');
      fprintf('输入 "nan" 表示未知或不适用。\n');

      % Read input line by line for upper triangle
      tempMatrix = NaN(n,n);
      try
          row_num = 1;
          while row_num < n
               prompt = sprintf('第 %d 行 (节点 %d 到 %d..%d): ', row_num, row_num, row_num+1, n);
               line_input = input(prompt, 's');
               if isempty(strtrim(line_input))
                   % Assume remaining elements for this row are NaN/Not specified
                   % We only need values up to column n anyway
               else
                  vals = str2num(line_input); %#ok<ST2NM> Use str2num cautiously
                  num_expected = n - row_num;
                  if length(vals) > num_expected
                       fprintf('错误：第 %d 行输入的值过多 (最多 %d 个)。\n', row_num, num_expected);
                       error('输入错误'); % Trigger catch block
                  end
                  % Place values into the correct columns (row_num+1 to n)
                  tempMatrix(row_num, row_num+1 : row_num+length(vals)) = vals;
               end
               row_num = row_num + 1;
          end

          % Validate input values based on type
          for r = 1:n
              for c = r+1:n
                  val = tempMatrix(r, c);
                  if strcmp(matrixName, 'Connectivity')
                      if ~isnan(val) && val ~= 0 && val ~= 1
                          fprintf('错误: 连接矩阵元素 (%d, %d) 必须是 0, 1, 或 NaN。\n', r, c);
                          error('输入错误');
                      elseif isnan(val) % Default unspecified connectivity to 0
                           tempMatrix(r,c) = 0;
                      end
                  elseif strcmp(matrixName, 'Resistance')
                       if ~isnan(val) && (val < allowedRange(1) || (~isinf(val) && val > allowedRange(2)))
                            fprintf('错误: 电阻值 (%d, %d) 必须 >= 0。\n', r, c);
                            error('输入错误');
                       end
                  else % Voltage or Current
                       if ~isnan(val) && (val < allowedRange(1) || val > allowedRange(2))
                            fprintf('错误: %s 值 (%d, %d) 超出范围。\n', matrixName, r, c);
                            error('输入错误');
                       end
                  end
              end
          end

          matrix = tempMatrix;
          validInput = true;

      catch ME
          fprintf('输入处理错误: %s\n', ME.message);
          fprintf('请确保使用数字、NaN、空格和分号正确输入。\n');
          % Ask user if they want to retry input for this matrix
          retry = input('是否重试输入此矩阵? (y/n): ', 's');
          if lower(retry) ~= 'y'
              matrix = []; % Indicate failure
              return; % Exit function, main loop will restart
          end
          % If retry is 'y', the loop continues
      end
  end
end

% --- 辅助函数：显示带有 'x' 的矩阵 ---
function displayMatrixWithX(matrix)
  [rows, cols] = size(matrix);
  strMatrix = cell(rows, cols); % Use cell array to store strings

  % Populate the cell array with formatted strings or 'x'
  for r = 1:rows
      for c = 1:cols
          if r >= c % Lower triangle and diagonal
              strMatrix{r, c} = 'x';
          else
              val = matrix(r, c);
              if isnan(val)
                  strMatrix{r, c} = 'NaN'; % Represent NaN
              elseif isinf(val)
                  % --- Corrected section for Infinity ---
                  if sign(val) == 1
                      signChar = '+';
                  else
                      signChar = '-';
                  end
                  strMatrix{r, c} = sprintf('%sInf', signChar); % Format as +Inf or -Inf
                  % --- End of correction ---
              else
                  strMatrix{r, c} = sprintf('%.4f', val); % Format regular numbers
              end
          end
      end
  end

  % Find the maximum width needed for each column for alignment
  maxWidths = zeros(1, cols);
  for c = 1:cols
      for r = 1:rows
          maxWidths(c) = max(maxWidths(c), length(strMatrix{r,c}));
      end
  end

  % Print the formatted matrix with alignment
  for r = 1:rows
      fprintf('  '); % Indent each row
      for c = 1:cols
          % Use dynamic field width for right alignment
          fprintf(['%' num2str(maxWidths(c)) 's  '], strMatrix{r,c});
      end
      fprintf('\n'); % New line after each row
  end
end

% --- 运行主函数 ---
% circuitSolverMNA(); % Call the function (if not running as a script)