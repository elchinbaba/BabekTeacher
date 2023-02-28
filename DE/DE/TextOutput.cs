namespace DE
{
	public class TextOutput
	{
		public delegate void WriteDelegate(object str);

		internal static void WriteLocal(object str)
		{
			System.Console.Write(str);
		}

		internal static void WriteLineLocal(object str)
		{
			System.Console.WriteLine(str);
		}

		public static WriteDelegate Write = WriteLocal;
		public static WriteDelegate WriteLine = WriteLineLocal;
	}
}
