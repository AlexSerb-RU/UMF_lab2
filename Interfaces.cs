using System;
using System.Collections.Generic;
using System.Text;

namespace UMF_lab2;

public interface IBasisFunction
{
    double Value(double x);
    double Derivative(double x);
}

public interface IFiniteElement
{
    int Id { get; }
    IReadOnlyList<int> GlobalNodes { get; }  // Глобальные индексы узлов (левый, правый)
    IReadOnlyList<IBasisFunction> BasisFunctions { get; }
    double Length { get; }
}

public interface IMesh
{
    IReadOnlyList<double> Nodes { get; }
    IReadOnlyList<IFiniteElement> Elements { get; }
    int NodeCount { get; }
    int ElementCount { get; }
    double GetNodeCoordinate(int nodeId);
    IFiniteElement GetElement(int elementId);
}

public interface ICoefficients
{


    double Lambda(double u, double x, double t = 0);

    double Sigma(double u, double x, double dudx, double t = 0);

    double SigmaDuDt(double u, double x, double t = 0);

    double F(double u, double x, double t = 0);

    // для краевого условия 2-го рода
    double Theta(double u, double x, double t = 0);
}

public interface IBoundaryConditions
{
    // Тип условий
    BoundaryType LeftType { get; }
    BoundaryType RightType { get; }

    // значения на границе (для 1-го рода)
    double LeftValue { get; }
    double RightValue { get; }

    // Функции для граничных условий
    Func<double, double, double> Theta { get; }

    bool IsBoundaryNode(int nodeId, int totalNodes);
}

public enum BoundaryType
{
    Dirichlet,   // 1-го рода
    Neumann,     // 2-го рода
    Robin        // 3-го рода
}

public interface ILocalAssembly
{
    double[,] LocalMatrix { get; }
    double[] LocalVector { get; }
    int Size { get; }

    void Compute(IFiniteElement element, ICoefficients coeffs,
                    IReadOnlyList<double> globalSolution);
}

public interface IGlobalAssembly
{
    BandedMatrix GlobalMatrix { get; }
    double[] GlobalVector { get; }

    void Assemble(IMesh mesh, ILocalAssembly localAssembly,
        ICoefficients coeffs, IReadOnlyList<double> solution,
        IBoundaryConditions bc);

    void ApplyBoundaryConditions(IBoundaryConditions bc, IReadOnlyList<double> solution);
}

public interface ILinearSolver
{
    double[] Solve(BandedMatrix matrix, double[] rhs);
}

public interface INonlinearSolver
{
    NonlinearSolverResult Solve(INonlinearProblem problem, double[] initialGuess);
}

public interface INonlinearProblem
{
    void ComputeSystem(IReadOnlyList<double> q, out BandedMatrix A, out double[] b);

    double[] ComputeResidual(IReadOnlyList<double> q);

    double ComputeResidualNorm(IReadOnlyList<double> q);

    int Size { get; }
}

public interface ILinearizationStrategy
{
    void Linearize(ILocalAssembly localAssembly, IFiniteElement element,
        ICoefficients coeffs, IReadOnlyList<double> solution,
        out double[,] linearizedMatrix, out double[] linearizedVector);
}

public interface IQuadrature
{
    IReadOnlyList<double> Points { get; }
    IReadOnlyList<double> Weights { get; }
    int Order { get; }
}

public interface IElementFactory
{
    IFiniteElement CreateElement(int id, int node1, int node2);
    IBasisFunction CreateLinearBasisFunction(int localId, double x1, double x2);
    IBasisFunction CreateQuadraticBasisFunction(int localId, double x1, double x2, double x3);
}

public interface ISolutionResult
{
    double[] Solution { get; }
    double[] ResidualHistory { get; }
    int IterationsCount { get; }
    bool Converged { get; }
    double FinalResidual { get; }
    double ComputationTime { get; }
}
