// DE Optimization Algo:rithm
// 15:11 4/15/2021

// Changes:
// 4:07 PM 6/30/2022: TextOutput class moved to DE namespace scope...
// 12:58 AM 7/2/2022: Addtional properties or Problem class: IsFeasible, Cost, Error
// 6:53 PM 7/4/2022: Problem.ResetSolution() added
// 5:42 PM 7/5/2022: Initialize() is added. Behaviour of Run() has changed
// 5:19 PM 7/14/2022: Problem.ResetSoltion(Ranges) and other overloads added
// 11:37 AM 11/1/2022: CostFun and ErrorFun marked protected

namespace DE
{
	public delegate double Function(double[] x);
	public delegate bool ExitFunction();
}
