using System;
using UnityEngine;

namespace PhysiK.Unity
{
    public sealed class PhysikTetMeshAuthoring : MonoBehaviour
    {
        public PhysikMaterialAsset material;

        private IntPtr world;
        private PhysikComponentHandle component;

        public PhysikComponentHandle CreateNativeComponent(
            IntPtr worldHandle,
            int[] nodeIndices,
            int[] tetNodeIndices)
        {
            if (worldHandle == IntPtr.Zero)
            {
                component = default;
                return component;
            }

            if (nodeIndices == null || tetNodeIndices == null || material == null)
            {
                component = default;
                return component;
            }

            world = worldHandle;
            PhysikMaterialDesc nativeMaterial = material.ToNative();
            component = PhysikNative.PHYSIK_CreateTetMeshComponentWithMaterialDesc(
                world,
                nodeIndices,
                nodeIndices.Length,
                tetNodeIndices,
                tetNodeIndices.Length / 4,
                ref nativeMaterial);
            return component;
        }

        public void ApplyMaterial()
        {
            if (world == IntPtr.Zero || !component.IsValid || material == null)
            {
                return;
            }

            PhysikMaterialDesc nativeMaterial = material.ToNative();
            PhysikNative.PHYSIK_SetTetMeshMaterial(
                world,
                component,
                ref nativeMaterial);
        }
    }
}
