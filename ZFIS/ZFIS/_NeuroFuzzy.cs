namespace ZFIS
{
	public class _NeuroFuzzy
	{
		int NRules;
		internal int NLParams;
		internal int NRParams;
		internal _Parameter[] LParams;
		internal _Parameter[] RParams;
		_Rule[] Rules;
		double[] Output;
		InferenceEngine This;

		public _NeuroFuzzy(_Parameter[] lps, _Parameter[] rps, _Rule[] rules)
		{
			NLParams = lps.Length;
			NRParams = rps.Length;
			NRules = rules.Length;
			LParams = lps;
			RParams = rps;
			Rules = rules;
			string errmsg;
			if ((errmsg = CheckIntegrity()) != NoErrors) System.Console.WriteLine(errmsg);
			Output = new double[NRParams];
			This = new InferenceEngine(this);
		}

		public const string NoErrors = "";

		string CheckIntegrity()
		{
			int maxn = 1;
			int i;
			string msg = "Integrity Checker: ";
			for (i = 0; i < NLParams; i++) maxn *= LParams[i].NTerms;
			if (NRules > maxn) return msg + "error: more than one rule with the same antecedent part detected";
			bool[][] termuse = new bool[NLParams][];
			for (i = 0; i < NLParams; i++) termuse[i] = new bool[LParams[i].NTerms];
			int r;
			int p;
			for (r = 0; r < NRules; r++) for (p = 0; p < NLParams; p++)
				{
					termuse[p][Rules[r].LeftSideTermIndexes[p]] = true;
				}
			for (p = 0; p < NLParams; p++) for (i = 0; i < termuse[p].Length; i++)
					if (!termuse[p][i]) return msg + "warning: some antecedent terms remain unused in rules";
			return NoErrors;
		}

		public double[] DoInference(double[] input)
		{
			double[] Output = new double[NRParams];
			This.Inference(input);
			int p;
			for (p = 0; p < NRParams; p++) Output[p] = this.Output[p];
			return Output;
		}

		// Inference Methods

		class InferenceEngine
		{
			// Takagi
			double[] outputL;
			double[] outputR;
			// Mamdani
			FuzzySet[] foutputL;
			FuzzySet[] foutputR;

			_NeuroFuzzy NF;

			internal delegate void InferenceRef(double[] input);
			internal InferenceRef Inference;

			internal InferenceEngine(_NeuroFuzzy nf)
			{
				NF = nf;
				outputL = new double[NF.NRParams];
				outputR = new double[NF.NRParams];
				foutputL = new FuzzySet[NF.NRParams];
				foutputR = new FuzzySet[NF.NRParams];
				Inference = Inference2;
			}

			//
			// Takagi approach
			void Inference1(double[] input)
			{
				int r;
				int p;
				double rf_sum_L = 0;
				double rf_sum_R = 0;
				for (p = 0; p < NF.NRParams; p++) outputL[p] = outputR[p] = 0;
				for (r = 0; r < NF.NRules; r++)
				{
					double rule_firing_L = 1;
					double rule_firing_R = 1;
					for (p = 0; p < NF.NLParams; p++)
					{
						rule_firing_L = System.Math.Min(
							rule_firing_L,
							NF.LParams[p].Terms[NF.Rules[r].LeftSideTermIndexes[p]].Value.Membership(input[p]).L
						);
						rule_firing_R = System.Math.Min(
							rule_firing_R,
							NF.LParams[p].Terms[NF.Rules[r].LeftSideTermIndexes[p]].Value.Membership(input[p]).R
						);
					}
					for (p = 0; p < NF.NRParams; p++)
					{
						outputL[p] +=
							rule_firing_L * NF.RParams[p].Terms[NF.Rules[r].RightSideTermIndexes[p]].Value.AL;
						outputR[p] +=
							rule_firing_R * NF.RParams[p].Terms[NF.Rules[r].RightSideTermIndexes[p]].Value.AR;
					}
					rf_sum_L += rule_firing_L;
					rf_sum_R += rule_firing_R;
				}
				if (rf_sum_L == 0) rf_sum_L = 1;
				if (rf_sum_R == 0) rf_sum_R = 1;
				for (p = 0; p < NF.NRParams; p++)
					NF.Output[p] = (outputL[p] / rf_sum_L + outputR[p] / rf_sum_R) / 2;
			}
			//

			// Mamdani approach
			void Inference2(double[] input)
			{
				int r;
				int p;
				for (p = 0; p < NF.NRParams; p++) foutputL[p] = foutputR[p] = FuzzySet.EmptySet;
				for (p = 0; p < NF.NRParams; p++) NF.Output[p] = 0;
				for (r = 0; r < NF.NRules; r++)
				{
					double rule_firing_L = 1;
					double rule_firing_R = 1;
					for (p = 0; p < NF.NLParams; p++)
					{
						rule_firing_L = System.Math.Min(
							rule_firing_L,
							NF.LParams[p].Terms[NF.Rules[r].LeftSideTermIndexes[p]].Value.Membership(input[p]).L
						);
						rule_firing_R = System.Math.Min(
							rule_firing_R,
							NF.LParams[p].Terms[NF.Rules[r].LeftSideTermIndexes[p]].Value.Membership(input[p]).R
						);
					}
					for (p = 0; p < NF.NRParams; p++)
					{
						FuzzySet s = NF.RParams[p].Terms[NF.Rules[r].RightSideTermIndexes[p]].Value.ToFuzzySet();
						foutputL[p] |= s.BoundMembershipToRange(0, rule_firing_L); // note: this method generates new capped fuzzy set w/o modifying s!
						foutputR[p] |= s.BoundMembershipToRange(0, rule_firing_R);
					}
				}
				for (p = 0; p < NF.NRParams; p++)
				{
					//Output[p]= (foutputL[p].Defuzzify()+foutputR[p].Defuzzify())/2;
					NF.Output[p] = (foutputL[p] | foutputR[p]).Defuzzify();
				}
			}
		}

		public double MaxError(double[][] inputs, double[][] desoutputs)
		{
			int no = desoutputs[0].Length;
			int ns = inputs.Length;
			//System.Console.WriteLine(no+" "+ns); System.Console.Read();
			int i, s;
			double err = 0.0;
			for (s = 0; s < ns; s++)
			{
				This.Inference(inputs[s]);
				for (i = 0; i < no; i++)
				{
					err = System.Math.Max(err, System.Math.Abs(Output[i] - desoutputs[s][i]));
				}
			}
			return err;
		}

		public double MSE(double[][] inputs, double[][] desoutputs)
		{
			int no = desoutputs[0].Length;
			int ns = inputs.Length;
			int i, s;
			double err = 0.0;
			for (s = 0; s < ns; s++)
			{
				This.Inference(inputs[s]);
				for (i = 0; i < no; i++)
				{
					double e = (Output[i] - desoutputs[s][i]);
					err += e * e;
				}
			}
			return err / (no * ns);
		}

		public double RMSE(double[][] inputs, double[][] desoutputs)
		{
			return System.Math.Sqrt(MSE(inputs, desoutputs));
		}

		public void DoTrain(double[][] inputs, double[][] desoutputs, double[,] lranges, double[,] rranges)
		{
			_NeuroFuzzyTrainDE det = new _NeuroFuzzyTrainDE(this, inputs, desoutputs, lranges, rranges);
			new DE.Optimizer(det).Run(
				1000,
				new DE.SearchControls(det.parameter_number * 100, 0.6, 0.9),
				new DE.BestStrategy(),
				new DE.Ranges(det.rmin, det.rmax)
			);
			det.LoadBest();
			det.Save("denfopt.txt");
		}

		public void DoTrain(double[][] inputs, double[][] desoutputs, double[,] lranges, double[,] rranges, int iterations, int strategy_no)
		{
			_NeuroFuzzyTrainDE det = new _NeuroFuzzyTrainDE(this, inputs, desoutputs, lranges, rranges);
			DE.Optimizer deo = new DE.Optimizer(det);
			DE.SearchStrategy strategy;
			if (strategy_no == 0) strategy = new DE.BestStrategy(); else strategy = new DE.RandStrategy();
			deo.Run(
				iterations,
				new DE.SearchControls(det.parameter_number * 10, 0.6, 0.9),
				strategy,
				new DE.Ranges(det.rmin, det.rmax)
			);
			det.LoadBest();
			det.Save("denfopt.txt");
		}

		override public string ToString()
		{
			string str = "";
			int r, p;
			for (r = 0; r < NRules; r++)
			{
				str += "IF ";
				for (p = 0; p < NLParams; p++)
				{
					str += LParams[p].Label + " IS " + "\"" + LParams[p].Terms[Rules[r].LeftSideTermIndexes[p]].Label + "\"";
					if (p < NLParams - 1) str += " AND ";
				}
				str += "\nTHEN ";
				for (p = 0; p < NRParams; p++)
				{
					//str+= RParams[p].Label+" IS "+RParams[p].Terms[Rules[r].RightSideTermIndexes[p]].Value;
					str += RParams[p].Label + " IS " + "\"" + RParams[p].Terms[Rules[r].RightSideTermIndexes[p]].Label + "\"";
					if (p < NRParams - 1) str += " AND ";
				}
				str += "\n";
			}
			return str;
		}

		public string ToString(double[] input)
		{
			if (input.Length < NLParams) return ToString();
			string str = "";
			int r, p;
			for (r = 0; r < NRules; r++)
			{
				str += "IF ";
				for (p = 0; p < NLParams; p++)
				{
					str += LParams[p].Label + " IS " + "\"" + LParams[p].Terms[Rules[r].LeftSideTermIndexes[p]].Label + "\"";
					str += "{" + LParams[p].Terms[Rules[r].LeftSideTermIndexes[p]].Value.Membership(input[p]) + "}";
					if (p < NLParams - 1) str += " AND ";
				}
				str += "\nTHEN ";
				for (p = 0; p < NRParams; p++)
				{
					//str+= RParams[p].Label+" IS "+RParams[p].Terms[Rules[r].RightSideTermIndexes[p]].Value;
					str += RParams[p].Label + " IS " + "\"" + RParams[p].Terms[Rules[r].RightSideTermIndexes[p]].Label + "\"";
					if (p < NRParams - 1) str += " AND ";
				}
				str += "\n";
			}
			This.Inference(input);
			str += "Final input: \n";
			for (p = 0; p < NLParams; p++)
			{
				str += LParams[p].Label + " = " + input[p] + "\n";
			}
			str += "Final Output: \n";
			for (p = 0; p < NRParams; p++)
			{
				str += RParams[p].Label + " = " + Output[p] + "\n";
			}
			return str;
		}
	}
}
