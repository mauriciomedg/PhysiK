using System;
using System.Runtime.InteropServices;

namespace PhysiK.Unity
{
    [StructLayout(LayoutKind.Sequential)]
    public struct PhysikMaterialDesc
    {
        public float density;
        public float youngModulus;
        public float poissonRatio;
        public float damping;
    }

    [StructLayout(LayoutKind.Sequential)]
    public struct PhysikComponentHandle
    {
        public uint index;
        public uint generation;

        public bool IsValid => index != 0xFFFFFFFFu && generation != 0u;
    }

    public static class PhysikNative
    {
        private const string DllName = "PhysiK";

        [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
        public static extern PhysikComponentHandle PHYSIK_CreateTetMeshComponentWithMaterialDesc(
            IntPtr world,
            int[] nodeIndices,
            int nodeCount,
            int[] tetNodeIndices,
            int tetCount,
            ref PhysikMaterialDesc material);

        [DllImport(DllName, CallingConvention = CallingConvention.Cdecl)]
        public static extern void PHYSIK_SetTetMeshMaterial(
            IntPtr world,
            PhysikComponentHandle component,
            ref PhysikMaterialDesc material);
    }
}
