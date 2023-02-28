namespace ZFIS
{
	public class FuzzySet
	{
		double[] X;
		double[] M;
		int N;
		bool Empty = false;

		public double[] Value { get { return X; } }
		public double[] Membership { get { return M; } }

		double MaxVal(double[] x)
		{
			double v = x[0];
			int i;
			for (i = 1; i < x.Length; i++) v = System.Math.Max(v, x[i]);
			return v;
		}

		public FuzzySet(double[] x, double[] m)
		{
			N = System.Math.Max(x.Length, m.Length);
			X = new double[N];
			M = new double[N];
			int i;
			for (i = 0; i < x.Length; i++) X[i] = x[i];
			for (; i < N; i++) X[i] = MaxVal(x) + i + 1 - N;
			for (i = 0; i < m.Length; i++) M[i] = m[i];
			for (; i < N; i++) M[i] = 0;
		}

		FuzzySet(double[] x, double[] m, int n)
		{
			N = System.Math.Min(x.Length, m.Length); // replace Max by Min: see above
			N = System.Math.Min(N, n);
			if (N == 0)
			{
				N = 1;
				X = new double[1];
				M = new double[1];
				X[0] = 0;
				M[0] = 0;
				Empty = true;
				return;
			}
			X = new double[N];
			M = new double[N];
			int i;
			for (i = 0; i < N; i++)
			{
				X[i] = x[i];
				M[i] = m[i];
			}
		}

		public FuzzySet(double x)
		{
			N = 1;
			X = new double[1];
			M = new double[1];
			X[0] = x;
			M[0] = 1;
		}

		public FuzzySet()
		{
			N = 1;
			X = new double[1];
			M = new double[1];
			X[0] = 0;
			M[0] = 0;
			Empty = true;
		}

		public static readonly FuzzySet EmptySet = new FuzzySet();

		public static FuzzySet operator |(FuzzySet a, FuzzySet b)
		{
			if (a.Empty) return b;
			if (b.Empty) return a;
			double[] tx = new double[a.N + b.N];
			double[] tm = new double[a.N + b.N];
			int i;
			for (i = 0; i < a.N; i++)
			{
				tx[i] = a.X[i];
				tm[i] = a.M[i];
			}
			int j;
			for (j = 0; j < b.N; j++)
			{
				tx[i + j] = b.X[j];
				tm[i + j] = b.M[j];
			}
			FuzzySet tf = new FuzzySet(tx, tm);
			return tf;
		}

		public static FuzzySet operator &(FuzzySet a, FuzzySet b)
		{
			if (a.Empty || b.Empty) return EmptySet;
			double[] tx = new double[a.N];
			double[] tm = new double[a.N];
			int i, j, k = 0;
			for (i = 0; i < a.N; i++) for (j = 0; j < b.N; j++)
				{
					if (a.X[i] == b.X[j])
					{
						tx[k] = a.X[i];
						tm[k] = System.Math.Min(a.M[i], b.M[j]);
						k++;
					}
				}
			return new FuzzySet(tx, tm, k);
		}

		public bool Equals(FuzzySet a)
		{
			if (this.Empty && a.Empty) return true;
			if (this.Empty || a.Empty) return false;
			if (N != a.N) return false;
			int i;
			for (i = 0; i < N; i++) if ((this.X[i] != a.X[i]) || (this.M[i] != a.M[i])) return false;
			return true;
		}

		public FuzzySet BoundMembershipToRange(double ml, double mr)
		{
			double[] tm = new double[N];
			for (int i = 0; i < N; i++)
			{
				tm[i] = System.Math.Max(System.Math.Min(M[i], mr), ml);
			}
			return new FuzzySet(X, tm);
		}

		public bool Contains(FuzzySet a)
		{
			if (a.Empty) return true;
			else if (this.Empty) return false;
			return (this | a).Equals(this);
		}

		public bool Contains(double a)
		{
			if (this.Empty) return false;
			int i;
			for (i = 0; i < N; i++) if (X[i] == a && M[i] == 1) return true;
			return false;
		}

		public double DefuzzifyCont() // Center of Gravity (COG) defuzzification
		{
			int i;
			double v1 = 0;
			double v2 = 0;
			if (N <= 1) return X[0];
			for (i = 0; i < N - 1; i++)
			{
				v2 += (M[i] + M[i + 1]) * (X[i + 1] - X[i]);
			}
			if (v2 == 0) return X[0];
			for (i = 0; i < N - 1; i++)
			{
				double v;
				v = 2 * (M[i] * X[i] + M[i + 1] * X[i + 1]) + M[i] * X[i + 1] + M[i + 1] * X[i];
				v *= (X[i + 1] - X[i]);
				v1 += v;
			}
			return v1 / v2 / 3.0;
		}

		public double Defuzzify() // Center of Gravity (COG) defuzzification
		{
			int i;
			double v1 = 0;
			double v2 = 0;
			if (N <= 1) return X[0];
			for (i = 0; i < N; i++)
			{
				v2 += M[i];
			}
			if (v2 == 0) return X[0];
			for (i = 0; i < N; i++)
			{
				v1 += X[i] * M[i];
			}
			return v1 / v2;
		}

		override public string ToString()
		{
			int i;
			string s = "{ ";
			if (!Empty) for (i = 0; i < N; i++)
				{
					s += X[i] + "/" + M[i];
					if (i < N - 1) s += ", ";
				}
			s += " }";
			return s;
		}
	}
}
