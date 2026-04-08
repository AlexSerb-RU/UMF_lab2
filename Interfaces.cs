using System;
using System.Collections.Generic;

namespace UMF_lab2;

public interface IBasisFunction
{
    double Value(double x);
    double Derivative(double x);
}

public interface IFiniteElement
{
    int Id { get; }
    IReadOnlyList<int> GlobalNodes { get; }
    IReadOnlyList<IBasisFunction> BasisFunctions { get; }
    double Length { get; }
}

public interface IMesh
{
    int NodeCount { get; }
}

public interface ISpaceMesh : IMesh
{
    IReadOnlyList<double> Nodes { get; }
    IReadOnlyList<IFiniteElement> Elements { get; }
    int ElementCount { get; }
    double LeftBoundary { get; }
    double RightBoundary { get; }
    double GetNodeCoordinate(int nodeId);
    IFiniteElement GetElement(int elementId);
}

public interface ITimeMesh : IMesh
{
    IReadOnlyList<double> TimePoints { get; }
    int TimeStepCount { get; }
    double TimeStart { get; }
    double TimeEnd { get; }
    double GetTimeStep(int stepIndex);
    double GetTimePoint(int stepIndex);
}

public interface ICoefficients
{
    double Lambda(double x, double t = 0.0);
    double Sigma(double dudx, double x, double t = 0.0);
    double DSigmaDDuDx(double dudx, double x, double t = 0.0);
    double F(double x, double t = 0.0);
}

public enum BoundaryType
{
    Dirichlet,
    Neumann
}

public interface IBoundaryConditions
{
    BoundaryType LeftType { get; }
    BoundaryType RightType { get; }
    double LeftValue(double t);
    double RightValue(double t);
}

public interface ILinearSolver
{
    double[] Solve(BandedMatrix matrix, double[] rhs);
}

public interface INonlinearProblem
{
    int Size { get; }
    void ComputeSystem(IReadOnlyList<double> q, out BandedMatrix a, out double[] b);
    //double[] ComputeResidual(IReadOnlyList<double> q);
    double ComputeRelativeResidual(IReadOnlyList<double> q);
}

public interface ILinearizableNonlinearProblem : INonlinearProblem
{
    void ComputeLinearizeMatrixAndResidual(IReadOnlyList<double> q, out BandedMatrix jacobian, out double[] residual);
}

public interface INonlinearSolver
{
    NonlinearSolverResult Solve(INonlinearProblem problem, double[] initialGuess);
}

public interface IInitialCondition
{
    double Value(double x);
}

public interface ITimeIntegrator
{
    double[] TimeStep(ILinearizableNonlinearProblem problem, double[] initialGuess);
    double[] LastResidualHistory { get; }
    int LastIterationsCount { get; }
}

public interface IExactSolution
{
    double Value(double x, double t);
}

public interface IManufacturedSolution
{
    double U(double x, double t);
    double Ux(double x, double t);
    double Uxx(double x, double t);
    double Ut(double x, double t);
}
